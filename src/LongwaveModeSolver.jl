"""
Main program/scratch file
"""
module LongwaveModeSolver

import Base: @_inline_meta, @propagate_inbounds, @_propagate_inbounds_meta
import Base: ==
using Dates
using LinearAlgebra

using StaticArrays
using StaticArrays: promote_tuple_eltype, convert_ntuple
using OrdinaryDiffEq, DiffEqCallbacks  # DiffEqBase,
using Parameters
using ProgressMeter
using BetterExp  # TEMP faster exp() until merged into Base
using JSON3, StructTypes
using Interpolations

using PolynomialRoots: roots!

using RootsAndPoles
using ModifiedHankelFunctionsOfOrderOneThird

# LongwaveModeSolver.jl
export bpm

# Geophysics.jl
export BField, Species, EigenAngle, Fields, Ground
export dip, azimuth
export waitprofile, electroncollisionfrequency, ioncollisionfrequency

# IO.jl
export BasicInput, TableInput, BatchInput, BasicOutput, BatchOutput

# Samplers.jl
export Receiver, GroundSampler

# Emitters.jl
export Transmitter, Dipole, VerticalDipole, HorizontalDipole, Frequency
export inclination, azimuth, altitude, fieldcomponent

# TMatrix.jl
export TMatrix

# Waveguides.jl
export HomogeneousWaveguide
export eigenangles

const BOTTOMHEIGHT = 0.0  # WARNING: if this isn't 0, many assumptions break

struct Dθ end

"""
    IntegrationParams{T}

Parameters passed to `OrdinaryDiffEq.jl` during the integration of the ionosphere reflection
coefficient matrix in `modefinder.jl`.

Fields:

    - tolerance::Float64
    - solver::T
    - force_dtmin::Bool

By default, the private function `integratedreflection` called by `findmodes` uses
`IntegrationParams(1e-6, OwrenZen5(), false)`.
"""
struct IntegrationParams{T}
    tolerance::Float64
    solver::T
    force_dtmin::Bool
end

export IntegrationParams

const DEFAULT_GRPFPARAMS = GRPFParams(100000, 1e-5, true)
const DEFAULT_INTEGRATIONPARAMS = IntegrationParams(1e-7, OwrenZen5(), false)

"""
    LWMSParams{T,H<:AbstractRange{Float64}}

Parameters for the `LongwaveModeSolver` module with defaults.

The struct is created with `Parameters.jl` `@with_kw` and supports that package's
instantiation capabilities, e.g.

```jldoctest
p = LWMSParams()
p2 = LWMSParams(earth_radius=6370e3)
p3 = LWMSParams(p2; grpf_params=GRPFParams(100000, 1e-6, true))
```
"""
@with_kw struct LWMSParams{T,H<:AbstractRange{Float64}}
    topheight::Float64 = 110e3
    earthradius::Float64 = 6369e3  # m
    earthcurvature::Bool = true
    curvatureheight::Float64 = 50e3  # m
    grpfparams::GRPFParams = DEFAULT_GRPFPARAMS
    integrationparams::IntegrationParams{T} = DEFAULT_INTEGRATIONPARAMS
    wavefieldheights::H = range(topheight, 0, length=513)
end
export LWMSParams

#
include("Antennas.jl")
include("EigenAngles.jl")
include("Emitters.jl")
include("Geophysics.jl")
include("Samplers.jl")
include("TMatrix.jl")
include("Waveguides.jl")
include("Wavefields.jl")

include("romberg.jl")

include("modeconversion.jl")
include("modefinder.jl")
include("magnetoionic.jl")
include("modesum.jl")

include("IO.jl")


function defaultcoordinates(frequency)
    # TODO: get a better idea of frequency transition
    if frequency > 15000
        Zb = deg2rad(complex(30.0, -10.0))
        Ze = deg2rad(complex(89.9, 0.0))
        d = deg2rad(60)
        Δr = deg2rad(0.4)
        coordgrid = eiwgdomain(Zb, Ze, d, Δr)
    else
        Zb = deg2rad(complex(0.0, -30.0))
        Ze = deg2rad(complex(89.9, 0.0))
        Δr = deg2rad(1.0)
        coordgrid = uppertriangularrectdomain(Zb, Ze, Δr)
    end

    return coordgrid
end
defaultcoordinates(f::Frequency) = defaultcoordinates(f.f)

"""
    bpm(waveguide::HomogeneousWaveguide, tx, rx::AbstractSampler{R}; coordgrid=nothing) where {R<:Real}

Specialized form for `AbstractSampler`s with a single distance which returns
scalar `E`, `amp`, and `phase`.
"""
function bpm(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler{R};
    coordgrid=nothing, params=LWMSParams()) where {R<:Real}

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

    E = Efield(modes, waveguide, tx, rx)
    amp = 10log10(abs2(E))  # == 20log10(abs(E))
    phase = angle(E)

    return E, amp, phase
end

"""
    bpm(waveguide::HomogeneousWaveguide, tx, rx; coordgrid=nothing)

Return electric field `E`, and field `amplitude` and `phase` using parameters:

    - `defaultcoordinates` for GRPF region
    - `PhysicalModeEquation`
"""
function bpm(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler;
    coordgrid=nothing, params=LWMSParams())

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

    E = Efield(modes, waveguide, tx, rx)

    amp = Vector{Float64}(undef, length(E))
    phase = similar(amp)
    @inbounds for i in eachindex(E)
        e = E[i]
        amp[i] = 10log10(abs2(e))  # == 20log10(abs(E))
        phase[i] = angle(e)  # ranges between -π:π rad
    end

    # Amplitude at transmitter may be calculated as Inf
    # TODO: replace with something accurate?
    isinf(amp[1]) && (amp[1] = 0)

    # By definition, phase at transmitter is 0, but is calculated as NaN
    isnan(phase[1]) && (phase[1] = 0)
    unwrap!(phase)

    return E, amp, phase
end


function bpm(waveguide::SegmentedWaveguide, tx::Emitter, rx::AbstractSampler;
    coordgrid=nothing, params=LWMSParams())

    if isnothing(coordgrid)
        coordgrid = defaultcoordinates(tx.frequency)
    else
        if minimum(imag(coordgrid)) < deg2rad(-31)
            @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields
                calculated with modified Hankel functions to overflow."
        end
    end

    numsegments = length(waveguide)

    heighttype = typeof(params.wavefieldheights)
    wavefields_vec = Vector{Wavefields{heighttype}}(undef, numsegments)
    adjwavefields_vec = Vector{Wavefields{heighttype}}(undef, numsegments)

    # Calculate wavefields and adjoint wavefields for each segment of waveguide
    for nsgmnt in 1:numsegments
        wvg = waveguide[nsgmnt]

        modeequation = PhysicalModeEquation(tx.frequency, wvg)

        modes = findmodes(modeequation, coordgrid, params=params)

        # adjoint wavefields are wavefields through adjoint waveguide, but for same modes
        # as wavefield
        @unpack bfield, species, ground = wvg
        adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
        adjwvg = HomogeneousWaveguide(adjoint_bfield, species, ground)

        wavefields = Wavefields(modes, params.wavefieldheights)
        adjwavefields = Wavefields(modes, params.wavefieldheights)

        calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg, params=params)

        wavefields_vec[nsgmnt] = wavefields
        adjwavefields_vec[nsgmnt] = adjwavefields
    end

    E = Efield(waveguide, wavefields_vec, adjwavefields_vec, tx, rx)

    amp = Vector{Float64}(undef, length(E))
    phase = similar(amp)
    @inbounds for i in eachindex(E)
        e = E[i]
        amp[i] = 10log10(abs2(e))  # == 20log10(abs(E))
        phase[i] = angle(e)  # ranges between -π:π rad
    end

    # Amplitude at transmitter may be calculated as Inf
    # TODO: replace with something accurate
    isinf(amp[1]) && (amp[1] = 0)

    # By definition, phase at transmitter is 0, but is calculated as NaN
    isnan(phase[1]) && (phase[1] = 0)
    unwrap!(phase)

    return E, amp, phase
end

"""
    bpm(filename::AbstractString; incrementalwrite=false, append=false, coordgrid=nothing)

Run the model given a String filename and save a JSON file of `BasicOutput`.
"""
function bpm(file::AbstractString;
    incrementalwrite=false, append=false,
    coordgrid=nothing, params=LWMSParams())

    ispath(file) || error("$file is not a valid file name")

    # Save output
    basepath = dirname(file)
    filename, fileextension = splitext(basename(file))

    outfile = joinpath(basepath, filename)*"_bpm.json"

    s = parse(file)
    if incrementalwrite
        s isa BatchInput || throw(ArgumentError("incrementalwrite only supported for BatchInput files"))
        output = buildrunsave(outfile, s,
                              append=append, coordgrid=coordgrid, params=params)
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
