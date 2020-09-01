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

"""
    BasicInput

# TODO: list default paramers that aren't specified in struct
"""
mutable struct BasicInput
    id::UUID
    name::String
    description::String
    datetime::DateTime

    # All units SI
    segment_ranges::Vector{Float64}
    hprimes::Vector{Float64}
    betas::Vector{Float64}
    b_mag::Vector{Float64}
    b_dip::Vector{Float64}
    b_az::Vector{Float64}
    ground_sigmas::Vector{Float64}
    ground_epsr::Vector{Int}
    frequency::Float64
    output_ranges::Vector{Float64}

    function BasicInput()
        s = new()
        setfield!(s, :frequency, NaN)
        return s
    end
end
StructTypes.StructType(::Type{BasicInput}) = StructTypes.Mutable()

mutable struct BasicOutput
    id::UUID
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
    iscomplete(s)

Return `true` if struct is completely defined, otherwise return `false`.
"""
function iscomplete(s)
    for fn in fieldnames(typeof(s))
        isdefined(s, fn) || return false
    end
    return true
end

"""
    validlengths(s)

Check if field lengths match.
"""
function validlengths(s)
    numsegments = length(s.segment_ranges)
    checkfields = (:hprimes, :betas, :b_mag, :b_dip, :b_az, :ground_sigmas, :ground_epsr)
    for field in checkfields
        length(getfield(s, field)) == numsegments || return false
    end
    return true
end

function parse(file)
    inputtypes = (BasicInput,)

    # To clarify the syntax here, `filecontents` is what is returned from inside
    # the `do` block; the json contents if it matched an input type, otherwise
    # nothing
    matched = false
    filecontents = open(file, "r") do f
        for t in inputtypes
            s = JSON3.read(f, t)
            if iscomplete(s) && validlengths(s)
                matched = true
                return s
            end
        end
        return nothing
    end

    # If we haven't returned, then no valid format
    if matched
        return filecontents
    else
        error("$file could not be matched to a valid format.")
    end
end

function buildandrun(s::BasicInput)
    if length(s.segment_ranges) == 1
        # HomogeneousWaveguide
        bfield = BField(only(s.b_mag), deg2rad(only(s.b_dip)), deg2rad(only(s.b_az)))
        species = Species(QE, ME, z -> waitprofile(z, only(s.hprimes), only(s.betas), cutoff_low=50e3, threshold=3e9),
                          electroncollisionfrequency)
        ground = Ground(only(s.ground_epsr), only(s.ground_sigmas))
        waveguide = HomogeneousWaveguide(bfield, species, ground)

        tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, FC_Ez)
    else
        # SegmentedWaveguide
        waveguide = SegmentedWaveguide(HomogeneousWaveguide)
        for i in eachindex(s.segment_ranges)
            bfield = BField(s.b_mag[i], deg2rad(s.b_dip[i]), deg2rad(s.b_az[i]))
            species = Species(QE, ME, z -> waitprofile(z, s.hprimes[i], s.betas[i], cutoff_low=50e3, threshold=3e9),
                              electroncollisionfrequency)
            ground = Ground(s.ground_epsr, s.ground_sigmas)
            push!(waveguide, HomogeneousWaveguide(bfield, species, ground, s.segment_ranges[i]))
        end
        tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, FC_Ez)
    end

    E, phase, amp = bpm(waveguide, tx, rx)

    return s.output_ranges, E, phase, amp
end
