"""
What information is needed?

Input:
    - N(z,r), ν(z,r), Bvec(r), σ(r), ϵᵣ(r)
    - In future, N(z,r,n) for n species
    - mass(n), charge(n) for n species

    - Should I interpolate N(z)?
    - Should I limit to h'(r), β(r)?
    - Should I build in conductivity map and allow defining places? (this is really an extra utility)

    - Additional parameters (earth radius, thresholds, etc)
        - Antenna power, orientation, frequency

    - Output info
        - What range(s)?
        - What E field component or amp, phase
        - (I need to fill in a receiver object)

Output:
    - Whatever is requested in input
        - Amp(r), phase(r)
        - E field component(r)?
"""

function jsonsafe!(v)
    for i in eachindex(v)
        if isnan(v[i]) || isinf(v[i])
            v[i] = 0
        end
    end
end

"""
    BasicInput
"""
mutable struct BasicInput
    name::String
    description::String
    datetime::DateTime

    # All units SI
    segment_ranges::Vector{Float64}
    hprimes::Vector{Float64}
    betas::Vector{Float64}
    b_mags::Vector{Float64}
    b_dips::Vector{Float64}
    b_azs::Vector{Float64}
    ground_sigmas::Vector{Float64}
    ground_epsrs::Vector{Int}
    frequency::Float64
    output_ranges::Vector{Float64}

    function BasicInput()
        s = new()
        setfield!(s, :frequency, NaN)
        return s
    end
end
StructTypes.StructType(::Type{BasicInput}) = StructTypes.Mutable()

mutable struct BatchInput{T}
    name::String
    description::String
    datetime::DateTime

    inputs::Vector{T}

    function BatchInput{T}() where T
        s = new{T}()
        return s
    end
end
BatchInput() = BatchInput{Any}()
StructTypes.StructType(::Type{<:BatchInput}) = StructTypes.Mutable()

mutable struct BasicOutput
    name::String
    description::String
    datetime::DateTime

    output_ranges::Vector{Float64}
    amplitude::Vector{Float64}
    phase::Vector{Float64}

    BasicOutput() = new()
end
StructTypes.StructType(::Type{BasicOutput}) = StructTypes.Mutable()

mutable struct BatchOutput{T}
    name::String
    description::String
    datetime::DateTime

    outputs::Vector{T}

    BatchOutput{T}() where {T} = new{T}()
end
BatchOutput() = BatchOutput{Any}()
StructTypes.StructType(::Type{<:BatchOutput}) = StructTypes.Mutable()
jsonsafe!(s::BatchOutput) = jsonsafe!(s.outputs)

"""
    iscomplete(s)

Return `true` if struct is completely defined, otherwise return `false`.
"""
function iscomplete(s)
    for fn in fieldnames(typeof(s))
        isdefined(s, fn) || return false
    end
    return true
end

function iscomplete(s::BatchInput)
    isdefined(s, :inputs) || return false
    for i in eachindex(s.inputs)
        iscomplete(s.inputs[i]) || return false
    end
    return true
end

"""
    validlengths(s)

Check if field lengths match.
"""
function validlengths(s)
    numsegments = length(s.segment_ranges)
    checkfields = (:hprimes, :betas, :b_mags, :b_dips, :b_azs, :ground_sigmas, :ground_epsrs)
    for field in checkfields
        length(getfield(s, field)) == numsegments || return false
    end
    return true
end

function validlengths(s::BatchInput)
    isdefined(s, :inputs) || return false
    for i in eachindex(s.inputs)
        validlengths(s.inputs[i]) || return false
    end
    return true
end

function parse(file)
    # More to less specific
    inputtypes = (BasicInput, BatchInput{BasicInput}, BatchInput{Any})

    matched = false
    let filecontents
        # Try supported types and accept it if it's valid
        for t in inputtypes
            # To clarify the syntax here, `filecontents` is what is returned from inside
            # the `do` block; the JSON contents or `nothing`
            filecontents = open(file, "r") do f
                # Once `f` is read, it can't be read again with another `t` so we open/close
                s = JSON3.read(f, t)
                if iscomplete(s) && validlengths(s)
                    matched = true
                    return s
                end
            end
            matched && break
        end

        if matched
            return filecontents
        else
            error("\"$file\" could not be matched to a valid format.")
        end
    end
end

function buildrun(s::BasicInput)
    if length(s.segment_ranges) == 1
        # HomogeneousWaveguide
        bfield = BField(only(s.b_mags), only(s.b_dips), only(s.b_azs))
        species = Species(QE, ME, z -> waitprofile(z, only(s.hprimes), only(s.betas), cutoff_low=50e3, threshold=3e9),
                          electroncollisionfrequency)
        ground = Ground(only(s.ground_epsrs), only(s.ground_sigmas))
        waveguide = HomogeneousWaveguide(bfield, species, ground)

        tx = Transmitter(VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, Fields.Ez)
    else
        # SegmentedWaveguide
        waveguide = SegmentedWaveguide(HomogeneousWaveguide)
        for i in eachindex(s.segment_ranges)
            bfield = BField(s.b_mags[i], s.b_dips[i], s.b_azs[i])
            species = Species(QE, ME, z -> waitprofile(z, s.hprimes[i], s.betas[i], cutoff_low=50e3, threshold=3e9),
                              electroncollisionfrequency)
            ground = Ground(s.ground_epsrs[i], s.ground_sigmas[i])
            push!(waveguide, HomogeneousWaveguide(bfield, species, ground, s.segment_ranges[i]))
        end
        tx = Transmitter(VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, Fields.Ez)
    end

    E, amp, phase = bpm(waveguide, tx, rx)

    output = BasicOutput()
    output.name = s.name
    output.description = s.description
    output.datetime = s.datetime

    output.output_ranges = s.output_ranges
    output.amplitude = amp
    output.phase = phase

    jsonsafe!(output.amp)
    jsonsafe!(output.phase)

    return output
end

function buildrun(s::BatchInput{BasicInput})
    batch = BatchOutput{BasicOutput}()
    batch.name = s.name
    batch.description = s.description
    batch.datetime = s.datetime

    for i in eachindex(s.inputs)
        output = buildandrun(s.inputs[i])
        push!(batch.outputs, output)
    end

    return batch
end

function buildrunsave(outfile, s::BatchInput{BasicInput}; append=false)
    if append
        batch = open(outfile, "r") do f
            s = JSON3.read(f, BatchOutput{BasicOutput})
            return s
        end
    else
        batch = BatchOutput{BasicOutput}()
        batch.name = s.name
        batch.description = s.description
        batch.datetime = s.datetime
    end

    for i in eachindex(s.inputs)
        name = s.inputs[i].name

        # Check if this case has already been run (useful for append)
        for o in eachindex(batch.outputs)
            name == batch.outputs[o].name && continue
        end

        output = buildandrun(s.inputs[i])
        push!(batch.outputs, output)

        json_str = JSON3.write(outputs)

        open(outfile*"_bpm.json", "w") do f
            write(f, json_str)
        end
    end

    return batch
end
