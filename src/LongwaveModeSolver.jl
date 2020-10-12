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
using BetterExp  # TEMP faster exp() until merged into Base
using JSON3, StructTypes
using Interpolations

using PolynomialRoots: roots!

using RootsAndPoles
using ModifiedHankelFunctionsOfOrderOneThird

using Debugger  # TEMP

# LongwaveModeSolver.jl
export bpm

# Geophysics.jl
export BField, Species, EigenAngle, Fields, Ground
export dip, azimuth
export waitprofile, electroncollisionfrequency, ioncollisionfrequency

# IO.jl
export BasicInput, BatchInput

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

#
const TOPHEIGHT = 110e3
const BOTTOMHEIGHT = zero(TOPHEIGHT)  # WARNING: if this isn't 0, many assumptions break

# Not great, but can be changed as `EARTHCURVATURE[]=false`
# TODO: where does this need to be considered?
const EARTHCURVATURE = Ref(true)

struct Dθ end

#
include("Antennas.jl")
include("EigenAngles.jl")
include("Emitters.jl")
include("Geophysics.jl")
include("IO.jl")
include("Samplers.jl")
include("TMatrix.jl")
include("Waveguides.jl")
include("Wavefields.jl")

include("romberg.jl")

include("magnetoionic.jl")
include("modeconversion.jl")
include("modefinder.jl")
include("modesum.jl")

function defaultcoordinates(frequency)
    # TODO: get a better idea of frequency transition
    if frequency > 15000
        Zb = deg2rad(complex(0.0, -10.0))
        Ze = deg2rad(complex(89.9, 0.0))
        d = deg2rad(60)
        Δr = deg2rad(0.5)
        origcoords = eiwgdomain(Zb, Ze, d, Δr)
    else
        Zb = deg2rad(complex(0.0, -30.0))
        Ze = deg2rad(complex(89.9, 0.0))
        Δr = deg2rad(1.0)
        origcoords = uppertriangularrectdomain(Zb, Ze, Δr)
    end

    return origcoords
end
defaultcoordinates(f::Frequency) = defaultcoordinates(f.f)

"""
    bpm(waveguide::HomogeneousWaveguide, tx, rx)

Return electric field `E`, and field `amplitude` and `phase` using parameters:

    - `defaultcoordinates` for GRPF region
    - `GRPFParams.tolerance = 1e-6`
    - `GRPFParams.multithreaded = true`
    - `susceptibilityinterpolator`
    - `PhysicalModeEquation`
"""
function bpm(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler)
    origcoords = defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = GRPFParams(est_num_nodes, 1e-6, true)

    Mfcn = susceptibilityinterpolator(tx.frequency, waveguide)
    modeequation = PhysicalModeEquation(tx.frequency, waveguide, Mfcn)

    E, amp, phase = bpm(waveguide, tx, rx, origcoords, grpfparams, modeequation)

    return E, amp, phase
end

"""
    bpm(waveguide::HomogeneousWaveguide, tx, rx::AbstractSampler{R}, origcoords, grpfparams,
    modeequation) where {R<:Real}

Specialized form for `AbstractSampler`s with a single distance which returns
scalar `E`, `amp`, and `phase`.
"""
function bpm(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler{R},
    origcoords, grpfparams::GRPFParams, modeequation::ModeEquation) where {R<:Real}
    if minimum(imag(origcoords)) < deg2rad(-31)
        @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields
            calculated with modified Hankel functions to overflow."
    end

    modes = findmodes(origcoords, grpfparams, modeequation)

    E = Efield(modes, waveguide, tx, rx)
    amp = 10log10(abs2(E))  # == 20log10(abs(E))
    phase = angle(E)

    return E, amp, phase
end

function bpm(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler,
    origcoords, grpfparams::GRPFParams, modeequation::ModeEquation)
    if minimum(imag(origcoords)) < deg2rad(-31)
        @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields
            calculated with modified Hankel functions to overflow."
    end

    modes = findmodes(origcoords, grpfparams, modeequation)

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

function bpm(waveguide::SegmentedWaveguide, tx, rx)
    origcoords = defaultcoordinates(tx.frequency)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = GRPFParams(est_num_nodes, 1e-6, true)

    E, amp, phase = bpm(waveguide, tx, rx, origcoords, grpfparams, true)
end

function bpm(waveguide::SegmentedWaveguide, tx::Emitter, rx::AbstractSampler,
    origcoords, grpfparams::GRPFParams, interpolateM::Bool)
    if minimum(imag(origcoords)) < deg2rad(-31)
        @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields
            calculated with modified Hankel functions to overflow."
    end


    numsegments = length(waveguide)

    # Predetermine types
    Mtype = eltype(susceptibility(TOPHEIGHT, tx.frequency, waveguide[1]))
    wftype = promote_type(Mtype, ComplexF64)  # ComplexF64 is type of EigenAngle

    wavefields_vec = Vector{Wavefields{wftype}}(undef, numsegments)
    adjwavefields_vec = Vector{Wavefields{wftype}}(undef, numsegments)

    # Calculate wavefields and adjoint wavefields for each segment of waveguide
    for nsgmnt in 1:numsegments
        wvg = waveguide[nsgmnt]

        if interpolateM
            # TODO: check if functionwrapper is necessary
            Mfcn = susceptibilityinterpolator(tx.frequency, wvg)
        else
            Mfcn = alt -> susceptibility(alt, tx.frequency, wvg)
        end
        modeequation = PhysicalModeEquation(tx.frequency, wvg, Mfcn)

        modes = findmodes(origcoords, grpfparams, modeequation)

        # adjoint wavefields are wavefields through adjoint waveguide, but for same modes
        # as wavefield
        @unpack bfield, species, ground = wvg
        adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
        adjwvg = HomogeneousWaveguide(adjoint_bfield, species, ground)

        # TODO< just empty and resize the Wavefields
        wavefields = Wavefields{wftype}(modes)
        adjwavefields = Wavefields{wftype}(modes)

        calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg)

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
    bpm(filename::AbstractString)

Run the model given a String filename and save a JSON file of `BasicOutput`.
"""
function bpm(file::AbstractString; incrementalwrite=false, append=false)
    ispath(file) || error("$file is not a valid file name")

    # Save output
    basepath = dirname(file)
    filename, fileextension = splitext(basename(file))

    outfile = joinpath(basepath, filename)

    s = parse(file)
    if incrementalwrite
        output = buildrunsave(outfile, s, append)
    else
        output = buildrun(s)

        json_str = JSON3.write(output)

        open(outfile*"_bpm.json", "w") do f
            write(f, json_str)
        end
    end

    return nothing
end

end
