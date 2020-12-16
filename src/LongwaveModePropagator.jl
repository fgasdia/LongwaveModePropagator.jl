module LongwaveModePropagator

import Base: @_inline_meta, @propagate_inbounds, @_propagate_inbounds_meta
import Base: ==
using Dates, LinearAlgebra

using StaticArrays
using StaticArrays: promote_tuple_eltype, convert_ntuple
using OrdinaryDiffEq, DiffEqCallbacks
using Parameters
using ProgressMeter
using BetterExp  # TEMP faster exp() until merged into Base
using JSON3, StructTypes
using Interpolations

using PolynomialRoots: roots!

using RootsAndPoles
using ModifiedHankelFunctionsOfOrderOneThird

# LongwaveModePropagator.jl
export propagate

# EigenAngle.jl
export attenuation, phasevelocity, referencetoground

# Geophysics.jl
export BField, Species, EigenAngle, Fields, Ground, GROUND
export waitprofile, electroncollisionfrequency, ioncollisionfrequency

# IO.jl
export BasicInput, TableInput, BatchInput, BasicOutput, BatchOutput

# Samplers.jl
export Receiver, Sampler, GroundSampler
export distance

# Emitters.jl
export Transmitter, Dipole, VerticalDipole, HorizontalDipole, Frequency
export inclination, azimuth, power

# modefinder.jl
export findmodes, PhysicalModeEquation

# Waveguides.jl
export HomogeneousWaveguide, SegmentedWaveguide


"""
    IntegrationParams{T}

Parameters passed to `OrdinaryDiffEq.jl` during the integration of the ionosphere reflection
coefficient matrix in `modefinder.jl`.

# Fields

- `tolerance::Float64 = 1e-7`: integration `atol` and `rtol`.
- `solver::T = OwrenZen5()`: a `DifferentialEquations.jl` solver.
- `force_dtmin::Bool = false`: if true, continue integration when solver reaches minimum
    step size.
"""
@with_kw struct IntegrationParams{T}
    tolerance::Float64 = 1e-7
    solver::T = OwrenZen5()
    force_dtmin::Bool = false
end
export IntegrationParams

const DEFAULT_GRPFPARAMS = GRPFParams(100000, 1e-5, true)

"""
    LMPParams{T,H <: AbstractRange{Float64}}

Parameters for the `LongwaveModePropagator` module with defaults:

- `topheight::Float64 = 110e3`: starting height for integration of the ionosphere reflection
    coefficient.
- `earthradius::Float64 = 6369`: Earth radius in meters.
- `earthcurvature::Bool = true`: toggle inclusion of Earth curvature in calculations. This
    is not supported by all functions.
- `curvatureheight::Float64 = 50e3`: reference height for Earth curvature in meters. At this
    height, the index of refraction is 1, and is therefore the reference height for
    eigenangles.
- `grpfparams::GRPFParams = GRPFParams(100000, 1e-5, true)`: parameters for the `GRPF`
    complex root-finding algorithm.
- `integrationparams::IntegrationParams{T} = IntegrationParams(1e-7, OwrenZen5(), false)`:
    parameters passed to `DifferentialEquations.jl` for integration of the ionosphere
    reflection coefficient.
- `wavefieldheights::H = range(topheight, 0, length=513)`: heights in meters at which
    wavefields will be integrated.

The struct is created using `Parameters.jl` `@with_kw` and supports that package's
instantiation capabilities, e.g.:

```jldoctest
p = LMPParams()
p2 = LMPParams(earth_radius=6370e3)
p3 = LMPParams(p2; grpf_params=GRPFParams(100000, 1e-6, true))
```
"""
@with_kw struct LMPParams{T,H<:AbstractRange{Float64}}
    topheight::Float64 = 110e3
    earthradius::Float64 = 6369e3  # m
    earthcurvature::Bool = true
    curvatureheight::Float64 = 50e3  # m
    grpfparams::GRPFParams = DEFAULT_GRPFPARAMS
    integrationparams::IntegrationParams{T} = IntegrationParams()
    wavefieldheights::H = range(topheight, 0, length=513)
end
export LMPParams

"""
    ModeEquation

Functions can dispatch on subtypes of `ModeEquation`, although currently only
`PhysicalModeEquation` is supported.

Future work might include the `ModifiedModeEquation` of [^Morfitt1976].

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA032573.
"""
abstract type ModeEquation end  # concrete subtypes defined in modefinder files

"Indicates derivative with respect to eigenangle ``θ``."
struct Dθ end

"""
Bottom of reflection coefficient integration. `LongwaveModePropagator.jl` assumes
`BOTTOMHEIGHT = 0.0`.
"""
const BOTTOMHEIGHT = 0.0  # WARNING: if this isn't zero, many assumptions break


#
include("utils.jl")
include("Antennas.jl")
include("Emitters.jl")  # must be before EigenAngles.jl
include("EigenAngles.jl")
include("Geophysics.jl")
include("Samplers.jl")
include("TMatrix.jl")
include("Waveguides.jl")
include("Wavefields.jl")

include("romberg.jl")

include("bookerquartic.jl")
include("magnetoionic.jl")
include("modeconversion.jl")
include("modefinder.jl")
include("modesum.jl")

include("IO.jl")


"""
    propagate(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler;
              coordgrid=nothing, params=LMPParams())

Compute electric field `E`, `amplitude`, and `phase`.
"""
function propagate(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler;
    coordgrid=nothing, params=LMPParams())

    if isnothing(coordgrid)
        coordgrid = defaultcoordinates(tx.frequency)
    else
        if minimum(imag(coordgrid)) < deg2rad(-31)
            @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields
                calculated with modified Hankel functions to overflow."
        end
    end

    modeequation = PhysicalModeEquation(tx.frequency, waveguide)
    modes = findmodes(modeequation, coordgrid, params=params)

    E = Efield(modes, waveguide, tx, rx, params=params)

    amplitude = Vector{Float64}(undef, length(E))
    phase = similar(amplitude)
    @inbounds for i in eachindex(E)
        e = E[i]
        amplitude[i] = 10log10(abs2(e))  # == 20log10(abs(E))
        phase[i] = angle(e)  # ranges between -π:π rad
    end

    unwrap!(phase)

    return E, amplitude, phase
end

"""
    propagate(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler{<:Real};
              coordgrid=nothing, params=LMPParams())

For `AbstractSampler`s with a single `distance`, compute scalar `E`, `amplitude`, and `phase`.
"""
function propagate(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler{<:Real};
    coordgrid=nothing, params=LMPParams())

    if isnothing(coordgrid)
        coordgrid = defaultcoordinates(tx.frequency)
    else
        if minimum(imag(coordgrid)) < deg2rad(-31)
            @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields
                calculated with modified Hankel functions to overflow."
        end
    end

    modeequation = PhysicalModeEquation(tx.frequency, waveguide)
    modes = findmodes(modeequation, coordgrid, params=params)

    E = Efield(modes, waveguide, tx, rx, params=params)
    amplitude = 10log10(abs2(E))  # == 20log10(abs(E))
    phase = angle(E)  # BUG? do we need to actually unwrap field all the way from tx?

    return E, amplitude, phase
end

"""
    propagate(waveguide::SegmentedWaveguide, tx::Emitter, rx::AbstractSampler;
              coordgrid=nothing, params=LMPParams())
"""
function propagate(waveguide::SegmentedWaveguide, tx::Emitter, rx::AbstractSampler;
    coordgrid=nothing, params=LMPParams())

    if isnothing(coordgrid)
        coordgrid = defaultcoordinates(tx.frequency)
    else
        if minimum(imag(coordgrid)) < deg2rad(-31)
            @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields
                calculated with modified Hankel functions to overflow."
        end
    end

    J = length(waveguide)  # number of waveguide segments

    heighttype = typeof(params.wavefieldheights)
    wavefields_vec = Vector{Wavefields{heighttype}}(undef, J)
    adjwavefields_vec = Vector{Wavefields{heighttype}}(undef, J)

    # Calculate wavefields and adjoint wavefields for each segment of waveguide
    for j in 1:J
        wvg = waveguide[j]

        modeequation = PhysicalModeEquation(tx.frequency, wvg)

        modes = findmodes(modeequation, coordgrid, params=params)

        # adjoint wavefields are wavefields through adjoint waveguide, but for same modes
        # as wavefield
        adjwvg = adjoint(wvg)

        wavefields = Wavefields(params.wavefieldheights, modes)
        adjwavefields = Wavefields(params.wavefieldheights, modes)

        calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg,
                              params=params)

        wavefields_vec[j] = wavefields
        adjwavefields_vec[j] = adjwavefields
    end

    E = Efield(waveguide, wavefields_vec, adjwavefields_vec, tx, rx, params=params)

    amplitude = Vector{Float64}(undef, length(E))
    phase = similar(amplitude)
    @inbounds for i in eachindex(E)
        e = E[i]
        amplitude[i] = 10log10(abs2(e))  # == 20log10(abs(E))
        phase[i] = angle(e)  # ranges between -π:π rad
    end

    unwrap!(phase)

    return E, amplitude, phase
end

"""
    propagate(file::AbstractString; outfile=missing, incrementalwrite=false, append=false,
              coordgrid=nothing)

Run the model scenario described by `file` and save the results as `outfile`.

If `outfile = missing`, the output file name will be `$(file)_output.json`.
"""
function propagate(file::AbstractString; outfile=missing, incrementalwrite=false,
    append=false, coordgrid=nothing, params=LMPParams())

    ispath(file) || error("$file is not a valid file name")

    if ismissing(outfile)
        basepath = dirname(file)
        filename, fileextension = splitext(basename(file))

        outfile = joinpath(basepath, filename)*"_output.json"
    end

    s = parse(file)
    if incrementalwrite
        s isa BatchInput || throw(ArgumentError("incrementalwrite only supported for BatchInput files"))
        output = buildrunsave(outfile, s, append=append, coordgrid=coordgrid, params=params)
    else
        output = buildrun(s, coordgrid=coordgrid, params=params)

        json_str = JSON3.write(output)

        open(outfile, "w") do f
            write(f, json_str)
        end
    end

    return output
end

end
