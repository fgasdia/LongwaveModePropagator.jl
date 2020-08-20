"""
Main program/scratch file
"""
module LongwaveModeSolver

using Printf  # TEMP

import Base: @_inline_meta, @propagate_inbounds, @_propagate_inbounds_meta
import Base: ==
using LinearAlgebra

using StaticArrays
using StaticArrays: promote_tuple_eltype, convert_ntuple
using DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks
using Parameters
# using NumericalIntegration
using BetterExp  # TEMP faster exp() until merged into Base

using PolynomialRoots: roots!

using RootsAndPoles
using ModifiedHankelFunctionsOfOrderOneThird

# Geophysics.jl
export BField, Species, EigenAngle, FieldComponent, Ground
export dip, azimuth
export waitprofile, electroncollisionfrequency, ioncollisionfrequency
export c₀, μ₀, ϵ₀, EARTH_RADIUS, CURVATURE_HEIGHT

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

#==
A struct like this could be used in place of the `const`s below.
It would allow us to parameterize the Complex type, but realistically it will
never be anything other than ComplexF64.

# PolynomialRoots package requires complex floats of arbitrary precision
struct BookerQuartic{T<:Complex{<:AbstractFloat}}
    roots::MVector{4,T}
    coeffs::MVector{5,T}
end
==#
# Passing MArrays between functions causes allocations. They are avoided by
# mutating this const in place. `roots!` requires Complex values and TMatrix
# is always ComplexF64, from which the coeffs are calculated.
const BOOKER_QUARTIC_ROOTS = MVector{4}(zeros(ComplexF64, 4))
const BOOKER_QUARTIC_COEFFS = MVector{5,ComplexF64}(undef)

# const BOOKER_QUARTIC_ROOTS = MVector{4}(zeros(Complex{BigFloat}, 4))
# const BOOKER_QUARTIC_COEFFS = MVector{5,Complex{BigFloat}}(undef)

struct Derivative_dθ end

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

    # By definition, phase at transmitter is 0, but is calculated as NaN
    isnan(phase[1]) && (phase[1] = 0)
    unwrap!(phase)

    return E, phase, amp
end

end
