"""
Main program/scratch file
"""
module LongwaveModeSolver

import Base: @_inline_meta, @propagate_inbounds, @_propagate_inbounds_meta
import Base: ==
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
export VACUUM_SPEED_OF_LIGHT, VACUUM_PERMEABILITY, VACUUM_PERMITTIVITY
export EARTH_RADIUS, CURVATURE_HEIGHT

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


"""
    triangulardomain(Za, Zb, Zc, Δr)

Generate initial mesh node coordinates for a *grid-aligned right* triangle domain
∈ {`Za`, `Zb`, `Zc`} with initial mesh step `Δr`.

This function generates a mesh grid for particular right triangles where two
sides of the triangle are aligned to the underlying real/imaginary grid. Examples
of such triangles are:

a ----- b
|     /
|   /
| /
c

where
    - a, b have greatest extent in x
    - a, c have greatest extent in y
"""
function triangulardomain(Za::Complex, Zb::Complex, Zc::Complex, Δr)
    rZa, iZa = reim(Za)
    rZb, iZb = reim(Zb)
    rZc, iZc = reim(Zc)

    #==
    # Check if this is a right triangle
    validtriangle = true
    if rZa == rZb == rZc
        validtriangle = false
    elseif iZa == iZb == iZc
        validtriangle = false
    elseif rZa == rZb
        iZa == iZc || iZb == iZc || (validtriangle = false)
    elseif rZa == rZc
        iZa == iZb || iZb == iZc || (validtriangle = false)
    elseif rZb == rZc
        iZa == iZb || iZa == iZc || (validtriangle = false)
    else
        validtriangle = false
    end
    validtriangle || throw(ArgumentError("`Za`, `Zb`, `Zc` do not define a grid-aligned right triangle"))

    iZa == iZb || ((Zb, Zc) = (Zc, Zb))
    rZb > rZa || ((Za, Zb) = (Zb, Za))
    ==#

    # Determine `dx` and `dy`
    X = rZb - rZa
    Y = abs(iZa - iZc)

    n = ceil(Int, Y/Δr + 1)
    dy = Y/(n-1)
    half_dy = dy/2

    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    # precalculate
    mn = m*n
    slope = Y/X

    vlength = mn + div(m, 2)  # BUG?
    T = promote_type(typeof(Za), typeof(Zb), typeof(Zc), typeof(Δr), Float64)
    v = Vector{T}()

    on = false
    for j = 0:m-1  # col
        for i = 0:n-1  # row (we're traversing down column)

            x = rZa + dx*j
            y = iZc + dy*i

            if (i+1) == n
                on = !on
            end
            if on
                y -= half_dy
            end

            if y >= (iZc + slope*(x - rZa))
                push!(v, complex(x, y))
            end
        end
    end

    return v
end

"""
    bpm(filename::AbstractString)

Run the model given a String filename.
"""
function bpm(file::AbstractString)
    ispath(file) || error("$filename is not a valid file name")

    s = parse(file)
    output_ranges, E, phase, amp = buildandrun(s)

    # Save output
    basepath = dirname(file)
    filename, fileextension = splitext(basename(file))

    output = BasicOutput()
    output.output_ranges = s.output_ranges
    output.amplitude = amp
    output.phase = phase

    json_str = JSON3.write(output)

    open(fullfile(basepath,filename*"_output",fileextension), "w") do f
        write(f, json_str)
    end
end

end
