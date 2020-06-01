using Test
using LinearAlgebra
using StaticArrays
using DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks
using RootsAndPoles

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C

# To profile in Juno
# @profiler (for i = 1:1000; LWMS.fcn(); end)

@testset "LongwaveModeSolver" begin
    include("test_geophysics.jl")
    include("test_wavefields.jl")
end
