using Test
using LinearAlgebra, Statistics
using StaticArrays
using DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks
using RootsAndPoles
using NLsolve
using JSON3, StructTypes

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C

# To profile in Juno
# @profiler (for i = 1:1000; LWMS.fcn(); end)

@testset "LongwaveModeSolver" begin
    include("test_EigenAngles.jl")
    include("test_Geophysics.jl")
    include("test_Waveguides.jl")

    include("test_wavefields.jl")
    include("test_modefinder.jl")

    include("test_IO.jl")
end
