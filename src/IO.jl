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

function parse(io)
    inputtypes = (BasicInput)
    open(io, "r") do f
        for t in inputtypes
            s = JSON3.read(f, t)
            completetype = iscomplete(s)
            if completetype
                validlengths(s) && return s
            end
        end
    end

    # If we haven't returned, then no valid format
    error("$io could not be matched to a valid format.")
end

function buildandrun(s::BasicInput)
    if length(s.segment_ranges) == 1
        # HomogeneousWaveguide
        bfield = BField(s.b_mag, deg2rad(s.b_dip), deg2rad(s.b_az))
        species = Species(qₑ, mₑ, z -> waitprofile(z, s.hprimes, s.betas),
                          electroncollisionfrequency)
        ground = Ground(s.epsr, s.sigmas)
        waveguide = HomogeneousWaveguide(bfield, species, ground)

        tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, LWMS.FC_Ez)
    else
        # SegmentedWaveguide
        waveguide = LWMS.SegmentedWaveguide(HomogeneousWaveguide)
        for i in eachindex(s.segment_ranges)
            bfield = BField(s.b_mag[i], deg2rad(s.b_dip[i]), deg2rad(s.b_az[i]))
            species = Species(qₑ, mₑ, z -> waitprofile(z, s.hprimes[i], s.betas[i]),
                              electroncollisionfrequency)
            ground = Ground(s.epsr, s.sigmas)
            push!(waveguide, HomogeneousWaveguide(bfield, species, ground))
        end
        tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(s.frequency), 100e3)
        rx = GroundSampler(s.output_ranges, LWMS.FC_Ez)
    end

    E, phase, amp = LWMS.bpm(waveguide, tx, rx)

    # TODO: we should be writing to file (if Julia input, no need to write to file)
    return s.output_ranges, E, phase, amp
end
