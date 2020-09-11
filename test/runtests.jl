using Test
using Dates
using LinearAlgebra, Statistics
using StaticArrays
using OrdinaryDiffEq, DiffEqCallbacks
using RootsAndPoles
using NLsolve
using JSON3, StructTypes

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const QE = LWMS.QE
const ME = LWMS.ME

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
