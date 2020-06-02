"""
Main program/scratch file
"""
module LongwaveModeSolver

import Base: @_inline_meta, @propagate_inbounds, @_propagate_inbounds_meta
using LinearAlgebra

using StaticArrays
using StaticArrays: promote_tuple_eltype, convert_ntuple
using DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks
using Parameters

using PolynomialRoots: roots!

using RootsAndPoles
using ModifiedHankelFunctionsOfOrderOneThird

# Geophysics.jl
export BField, Species, EigenAngle, FieldComponent, Ground
export dip, azimuth
export waitprofile, electroncollisionfrequency, ioncollisionfrequency
export Rₑ, c₀, μ₀, ϵ₀, CURVATURE_HEIGHT

# Scenario.jl
export Receiver, GroundSampler
export Transmitter, Dipole, VerticalDipole, HorizontalDipole, Frequency

# TMatrix.jl
export TMatrix

# PropagationModel.jl

#
const TOPHEIGHT = 95e3  # TODO: temporary - should be part of an actual IntegrationParameters
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
# mutating this const in place. `roots!` requires Complex values.
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

include("magnetoionic.jl")
include("modefinder.jl")
include("modesum.jl")
include("wavefields.jl")

end
