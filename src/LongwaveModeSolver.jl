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

using GRPF
using ModifiedHankelFunctionsOfOrderOneThird

# Geophysics.jl
export BField, Constituent, EigenAngle, FieldComponent, Ground
export waitprofile, electroncollisionfrequency, ioncollisionfrequency
export Rₑ, c₀, μ₀, ϵ₀, H

# Scenario.jl
export Receiver, GroundSampler
export Transmitter, Dipole, VerticalDipole, HorizontalDipole, Frequency

# TMatrix.jl
export TMatrix

# PropagationModel.jl


#
include("Geophysics.jl")
include("Scenario.jl")
include("TMatrix.jl")
include("PropagationModel.jl")

end
