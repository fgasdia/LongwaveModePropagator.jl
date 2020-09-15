"""
Main program/scratch file
"""
module LongwaveModeSolver

import Base: @_inline_meta, @propagate_inbounds, @_propagate_inbounds_meta
import Base: ==
using UUIDs, Dates
using LinearAlgebra

using StaticArrays
using StaticArrays: promote_tuple_eltype, convert_ntuple
using OrdinaryDiffEq, DiffEqCallbacks  # DiffEqBase,
using Parameters
using BetterExp  # TEMP faster exp() until merged into Base
using JSON3, StructTypes

using PolynomialRoots: roots!

using RootsAndPoles
using ModifiedHankelFunctionsOfOrderOneThird

# LongwaveModeSolver.jl
export bpm

# Geophysics.jl
export BField, Species, EigenAngle, FieldComponent, Ground
export dip, azimuth
export waitprofile, electroncollisionfrequency, ioncollisionfrequency

# IO.jl
export BasicInput

# Samplers.jl
export Receiver, GroundSampler

# Emitters.jl
export Transmitter, Dipole, VerticalDipole, HorizontalDipole, Frequency
export distance, elevation, azimuth, altitude, fieldcomponent

# TMatrix.jl
export TMatrix

# Waveguides.jl
export HomogeneousWaveguide
export eigenangles

#
const TOPHEIGHT = 110e3  # TODO: temporary - should be part of an actual IntegrationParameters
const BOTTOMHEIGHT = zero(TOPHEIGHT)  # should this be an actual const? Nothing works if it's not 0...

# Not great, but can be changed as `EARTHCURVATURE[]=false`
# TODO: where does this need to be considered?
const EARTHCURVATURE = Ref(true)

struct Derivative_dθ end

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


function bpm(waveguide::HomogeneousWaveguide, tx, rx)
    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8

    modes = findmodes(origcoords, tx.frequency, waveguide, tolerance)

    E = Efield(modes, waveguide, tx, rx)

    phase = Vector{Float64}(undef, length(E))
    amp = similar(phase)
    @inbounds for i in eachindex(E)
        e = E[i]
        phase[i] = angle(e)  # ranges between -π:π rad
        amp[i] = 10log10(abs2(e))  # == 20log10(abs(E))
    end

    # Amplitude at transmitter may be calculated as Inf
    # TODO: replace with something accurate?
    isinf(amp[1]) && (amp[1] = 0)

    # By definition, phase at transmitter is 0, but is calculated as NaN
    isnan(phase[1]) && (phase[1] = 0)
    unwrap!(phase)

    return E, phase, amp
end

function bpm(waveguide::SegmentedWaveguide, tx, rx)
    zs = range(TOPHEIGHT, 0, length=513)
    # zs = range(TOPHEIGHT, 0, length=129)
    nrsgmnt = length(waveguide)

    wavefields_vec = Vector{Wavefields{typeof(zs)}}(undef, nrsgmnt)
    adjwavefields_vec = Vector{Wavefields{typeof(zs)}}(undef, nrsgmnt)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-7

    for nsgmnt in 1:nrsgmnt
        wvg = waveguide[nsgmnt]
        modes = findmodes(origcoords, tx.frequency, wvg, tolerance)

        # TEMP: not necessary to sort, but easier to compare to LWPC
        sort!(modes,rev=true)
        # if nsgmnt > 1
        #     modes = modes[1:12]
        # end

        # adjoint wavefields are wavefields through adjoint waveguide, but for same modes as wavefield
        @unpack bfield, species, ground = wvg
        adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
        adjwvg = HomogeneousWaveguide(adjoint_bfield, species, ground)

        wavefields = Wavefields(modes, zs)
        adjwavefields = Wavefields(modes, zs)

        calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg)

        wavefields_vec[nsgmnt] = wavefields
        adjwavefields_vec[nsgmnt] = adjwavefields
    end

    E = Efield(waveguide, wavefields_vec, adjwavefields_vec, tx, rx)

    phase = Vector{Float64}(undef, length(E))
    amp = similar(phase)
    @inbounds for i in eachindex(E)
        e = E[i]
        phase[i] = angle(e)  # ranges between -π:π rad
        amp[i] = 10log10(abs2(e))  # == 20log10(abs(E))
    end

    # Amplitude at transmitter may be calculated as Inf
    # TODO: replace with something accurate?
    isinf(amp[1]) && (amp[1] = 0)

    # By definition, phase at transmitter is 0, but is calculated as NaN
    isnan(phase[1]) && (phase[1] = 0)
    unwrap!(phase)

    return E, phase, amp
end

"""
    bpm(filename::AbstractString)

Run the model given a String filename and save a JSON file of `BasicOutput`.
"""
function bpm(file::AbstractString)
    ispath(file) || error("$filename is not a valid file name")

    s = parse(file)
    output_ranges, E, phase, amp = buildandrun(s)

    # Save output
    basepath = dirname(file)
    filename, fileextension = splitext(basename(file))

    output = BasicOutput()
    output.name = s.name
    output.description = s.description
    output.datetime = s.datetime

    # TEMP, otherwise amp may be NaN which cannot be written to JSON
    amp[1] = amp[2]

    output.output_ranges = s.output_ranges
    output.amplitude = amp
    output.phase = rad2deg.(phase)

    json_str = JSON3.write(output)

    open(joinpath(basepath,filename*"_bpm"*fileextension), "w") do f
        write(f, json_str)
    end

    return nothing
end

end
