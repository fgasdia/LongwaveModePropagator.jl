#==
Functions related to running and saving LMP through JSON and other files.
==#

function jsonsafe!(v)
    for i in eachindex(v)
        if isnan(v[i]) || isinf(v[i])
            v[i] = 0
        end
    end
end

abstract type Input end

"""
    BasicInput

# Fields

- `name::String`
- `description::String`
- `datetime::DateTime`
- `segment_ranges::Vector{Float64}`: distance from transmitter to the beginning of each
    `HomogeneousWaveguide` segment in meters.
- `hprimes::Vector{Float64}`: Wait's ``h′`` parameter for each `HomogeneousWaveguide` segment.
- `betas::Vector{Float64}`: Wait's ``β`` parameter for each `HomogeneousWaveguide` segment.
- `b_mags::Vector{Float64}`: magnetic field magnitude for each `HomogeneousWaveguide` segment.
- `b_dips::Vector{Float64}`: magnetic field dip angles in radians for each
    `HomogeneousWaveguide` segment.
- `b_azs::Vector{Float64}`: magnetic field azimuth in radians "east" of the propagation
    direction for each `HomogeneousWaveguide` segment.
- `ground_sigmas::Vector{Float64}`: ground conductivity in Siemens per meter for each
    `HomogeneousWaveguide` segment.
- `ground_epsrs::Vector{Int}`: ground relative permittivity for each `HomogeneousWaveguide`
    segment.
- `frequency::Float64`: transmitter frequency in Hertz.
- `output_ranges::Vector{Float64}`: distances from the transmitter at which the field will
    be calculated.
"""
mutable struct BasicInput <: Input
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

"""
    TableInput <: Input

# Fields

- `name::String`
- `description::String`
- `datetime::DateTime`
- `segment_ranges::Vector{Float64}`: distance from transmitter to the beginning of each
    `HomogeneousWaveguide` segment in meters.
- `altitude::Vector{Float64}`: altitude above ground in meters for which the `density` and
    `collision_frequency` profiles are specified.
- `density::Vector{Float64}`: electron density at each `altitude` in ``m⁻³``.
- `collision_frequency::Vector{Float64}`: electron-ion collision frequency at each
    `altitude` in ``s⁻¹``.
- `b_dips::Vector{Float64}`: magnetic field dip angles in radians for each
    `HomogeneousWaveguide` segment.
- `b_azs::Vector{Float64}`: magnetic field azimuth in radians "east" of the propagation
    direction for each `HomogeneousWaveguide` segment.
- `ground_sigmas::Vector{Float64}`: ground conductivity in Siemens per meter for each
    `HomogeneousWaveguide` segment.
- `ground_epsrs::Vector{Int}`: ground relative permittivity for each `HomogeneousWaveguide`
    segment.
- `frequency::Float64`: transmitter frequency in Hertz.
- `output_ranges::Vector{Float64}`: distances from the transmitter at which the field will
    be calculated.
"""
mutable struct TableInput <: Input
    name::String
    description::String
    datetime::DateTime

    # All units SI
    segment_ranges::Vector{Float64}
    altitude::Vector{Float64}
    density::Vector{Vector{Float64}}
    collision_frequency::Vector{Vector{Float64}}
    b_mags::Vector{Float64}
    b_dips::Vector{Float64}
    b_azs::Vector{Float64}
    ground_sigmas::Vector{Float64}
    ground_epsrs::Vector{Int}
    frequency::Float64
    output_ranges::Vector{Float64}

    function TableInput()
        s = new()
        setfield!(s, :frequency, NaN)
        setfield!(s, :density, [Vector{Float64}()])
        setfield!(s, :collision_frequency, [Vector{Float64}()])
        return s
    end
end
StructTypes.StructType(::Type{TableInput}) = StructTypes.Mutable()

"""
    BatchInput{T} <: Input

A collection of `inputs` with a batch `name`, `description`, and `datetime`.
"""
mutable struct BatchInput{T} <: Input
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

abstract type Output end

"""
    BasicOutput <: Output

# Fields

- `name::String`
- `description::String`
- `datetime::DateTime`
- `output_ranges::Vector{Float64}`
- `amplitude::Vector{Float64}`
- `phase::Vector{Float64}`
"""
mutable struct BasicOutput <: Output
    name::String
    description::String
    datetime::DateTime

    output_ranges::Vector{Float64}
    amplitude::Vector{Float64}
    phase::Vector{Float64}

    BasicOutput() = new()
end
StructTypes.StructType(::Type{BasicOutput}) = StructTypes.Mutable()

"""
    BatchOutput{T} <: Output

A collection of `outputs` with a batch `name`, `description`, and `datetime`.

See also: [`BatchInput`](@ref)
"""
mutable struct BatchOutput{T} <: Output
    name::String
    description::String
    datetime::DateTime

    outputs::Vector{T}

    function BatchOutput{T}() where {T}
        s = new{T}()
        s.outputs = T[]
        return s
    end
end
BatchOutput() = BatchOutput{Any}()
StructTypes.StructType(::Type{<:BatchOutput}) = StructTypes.Mutable()
jsonsafe!(s::BatchOutput) = jsonsafe!(s.outputs)

"""
    iscomplete(s)

Return `true` if input or output struct `s` is completely defined, otherwise return `false`.
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

function iscomplete(s::BatchOutput)
    isdefined(s, :outputs) || return false
    for i in eachindex(s.outputs)
        iscomplete(s.outputs[i]) || return false
    end
    return true
end

"""
    validlengths(s)

Check if field lengths of input `s` match their number of segments.
"""
validlengths

function validlengths(s::BasicInput)
    numsegments = length(s.segment_ranges)
    checkfields = (:hprimes, :betas, :b_mags, :b_dips, :b_azs, :ground_sigmas,
        :ground_epsrs)
    for field in checkfields
        length(getfield(s, field)) == numsegments || return false
    end
    return true
end

function validlengths(s::TableInput)
    numsegments = length(s.segment_ranges)
    checkfields = (:b_mags, :b_dips, :b_azs, :ground_sigmas, :ground_epsrs)
    for field in checkfields
        length(getfield(s, field)) == numsegments || return false
    end

    numaltitudes = length(s.altitude)
    matrixfields = (:density, :collision_frequency)
    for field in matrixfields
        v = getfield(s, field)
        length(v) == numsegments || return false
        for i = 1:numsegments
            length(v[i]) == numaltitudes || return false
        end
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

"""
    parse(file)

Parse a JSON file compatible with `Input` or `Output` types.
"""
function parse(file)
    # More to less specific
    types = (BasicInput, TableInput,
        BatchInput{BasicInput}, BatchInput{TableInput}, BatchInput{Any},
        BasicOutput, BatchOutput{BasicOutput}, BatchOutput{Any})

    matched = false
    let filecontents
        for t in types
            filecontents = parse(file, t)
            if !isnothing(filecontents)
                matched = true
                break
            end
        end

        if matched
            return filecontents
        else
            error("\"$file\" could not be matched to a valid format.")
        end
    end
end

function parse(file, t::Type{<:Input})
    matched = false

    # To clarify the syntax here, `filecontents` is what is returned from inside
    # the `do` block; the JSON contents or `nothing`
    filecontents = open(file, "r") do f
        s = JSON3.read(f, t)
        if iscomplete(s) && validlengths(s)
            matched = true
            return s
        end
    end

    matched ? filecontents : nothing
end

function parse(file, t::Type{<:Output})
    matched = false

    # To clarify the syntax here, `filecontents` is what is returned from inside
    # the `do` block; the JSON contents or `nothing`
    filecontents = open(file, "r") do f
        s = JSON3.read(f, t)
        if iscomplete(s)
            matched = true
            return s
        end
    end

    matched ? filecontents : nothing
end

"""
    buildwaveguide(s::BasicInput, i)

Return `HomogeneousWaveguide` from the `i`th entry in each field of `s`.
"""
function buildwaveguide(s::BasicInput, i)
    bfield = BField(s.b_mags[i], s.b_dips[i], s.b_azs[i])
    species = Species(QE, ME, z -> waitprofile(z, s.hprimes[i], s.betas[i],
                                               cutoff_low=40e3),
                      electroncollisionfrequency)
    ground = Ground(s.ground_epsrs[i], s.ground_sigmas[i])
    return HomogeneousWaveguide(bfield, species, ground, s.segment_ranges[i])
end

"""
    buildwaveguide(s::TableInput, i)

Return `HomogeneousWaveguide` from the `i`th entry in each field of `s` with a linear
interpolation over `density` and `collision_frequency`.
"""
function buildwaveguide(s::TableInput, i)
    bfield = BField(s.b_mags[i], s.b_dips[i], s.b_azs[i])
    density_itp = LinearInterpolation(s.altitude, s.density[i], extrapolation_bc=Line())
    collision_itp = LinearInterpolation(s.altitude, s.collision_frequency[i], extrapolation_bc=Line())
    species = Species(QE, ME, density_itp, collision_itp)
    ground = Ground(s.ground_epsrs[i], s.ground_sigmas[i])
    return HomogeneousWaveguide(bfield, species, ground, s.segment_ranges[i])
end

"""
    buildrun(s::BasicInput; coordgrid=nothing, params=LMPParams())
    buildrun(s::TableInput; coordgrid=nothing, params=LMPParams())
    buildrun(s::BatchInput; coordgrid=nothing, params=LMPParams())

Build LMP structs from an `Input` and run `LMP`.
"""
buildrun

function buildrun(s::BasicInput; coordgrid=nothing, params=LMPParams())
    if length(s.segment_ranges) == 1
        # HomogeneousWaveguide
        bfield = BField(only(s.b_mags), only(s.b_dips), only(s.b_azs))
        species = Species(QE, ME, z -> waitprofile(z, only(s.hprimes), only(s.betas),
                                                   cutoff_low=40e3),
                          electroncollisionfrequency)
        ground = Ground(only(s.ground_epsrs), only(s.ground_sigmas))
        waveguide = HomogeneousWaveguide(bfield, species, ground)

        tx = Transmitter(VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, Fields.Ez)
    else
        # SegmentedWaveguide
        waveguide = SegmentedWaveguide([buildwaveguide(s, i) for i in
                                        eachindex(s.segment_ranges)])
        tx = Transmitter(VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, Fields.Ez)
    end

    E, amp, phase = propagate(waveguide, tx, rx, coordgrid=coordgrid, params=params)

    output = BasicOutput()
    output.name = s.name
    output.description = s.description
    output.datetime = Dates.now()

    output.output_ranges = s.output_ranges
    output.amplitude = amp
    output.phase = phase

    jsonsafe!(output.amplitude)
    jsonsafe!(output.phase)

    return output
end

function buildrun(s::TableInput; coordgrid=nothing, params=LMPParams())

    if length(s.segment_ranges) == 1
        # HomogeneousWaveguide
        bfield = BField(only(s.b_mags), only(s.b_dips), only(s.b_azs))
        density_itp = LinearInterpolation(s.altitude, only(s.density), extrapolation_bc=Line())
        collision_itp = LinearInterpolation(s.altitude, only(s.collision_frequency), extrapolation_bc=Line())
        species = Species(QE, ME, density_itp, collision_itp)
        ground = Ground(only(s.ground_epsrs), only(s.ground_sigmas))
        waveguide = HomogeneousWaveguide(bfield, species, ground)

        tx = Transmitter(VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, Fields.Ez)
    else
        # SegmentedWaveguide
        waveguide = SegmentedWaveguide([buildwaveguide(s, i) for i in
                                        eachindex(s.segment_ranges)])
        tx = Transmitter(VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, Fields.Ez)
    end

    E, amp, phase = propagate(waveguide, tx, rx, coordgrid=coordgrid, params=params)

    output = BasicOutput()
    output.name = s.name
    output.description = s.description
    output.datetime = Dates.now()

    output.output_ranges = s.output_ranges
    output.amplitude = amp
    output.phase = phase

    jsonsafe!(output.amplitude)
    jsonsafe!(output.phase)

    return output
end

function buildrun(s::BatchInput; coordgrid=nothing, params=LMPParams())

    batch = BatchOutput{BasicOutput}()
    batch.name = s.name
    batch.description = s.description
    batch.datetime = Dates.now()

    for i in eachindex(s.inputs)
        output = buildrun(s.inputs[i], coordgrid=coordgrid, params=params)
        push!(batch.outputs, output)
    end

    return batch
end

"""
    buildrunsave(outfile, s::BatchInput; append=false, coordgrid=nothing, params=LMPParams())

Similar to `buildrun`, except it saves results into `outfile` as `s` is processed.

If `append=true`, this function parses `outfile` for preexisting results and only runs the
remaining scenarios in `s`. Otherwise, a new `BatchOutput` is created.
"""
function buildrunsave(outfile, s::BatchInput; append=false, coordgrid=nothing,
    params=LMPParams())

    if append && isfile(outfile)
        batch = open(outfile, "r") do f
            v = JSON3.read(f, BatchOutput{BasicOutput})
            return v
        end
    else
        batch = BatchOutput{BasicOutput}()
        batch.name = s.name
        batch.description = s.description
        batch.datetime = Dates.now()
    end

    skip = false
    p = Progress(length(s.inputs), 5)
    for i in eachindex(s.inputs)
        name = s.inputs[i].name

        # Check if this case has already been run (useful for append)
        for o in eachindex(batch.outputs)
            if name == batch.outputs[o].name
                skip = true
                break
            end
        end
        if skip
            skip = false
            next!(p)
            continue
        end

        output = buildrun(s.inputs[i], coordgrid=coordgrid, params=params)
        push!(batch.outputs, output)

        json_str = JSON3.write(batch)

        open(outfile, "w") do f
            write(f, json_str)
        end

        next!(p)
    end

    return batch
end
