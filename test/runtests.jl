using Test
using Dates
using LinearAlgebra, Statistics
using StaticArrays, Parameters
using OrdinaryDiffEq, DiffEqCallbacks
using RootsAndPoles
using VoronoiDelaunay
using JSON3, StructTypes
using Interpolations, NLsolve, FiniteDiff

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const QE = LWMS.QE
const ME = LWMS.ME

# To profile in Juno
# @profiler (for i = 1:1000; LWMS.fcn(); end)

struct TestScenario{T1,T2,T3}
    ea::T1
    bfield::BField
    species::T2
    ground::Ground
    tx::T3
end

const verticalB_scenario = TestScenario(
    EigenAngle(deg2rad(complex(85.0, -1.0))),
    BField(50e-6, π/2, 0),
    Species(QE, ME,
            z->waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
            z->electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT)),
    Ground(15, 0.001),
    Transmitter(24e3))

const scenario = TestScenario(
    EigenAngle(1.4152764714690873 - 0.01755942938376613im),
    BField(50e-6, deg2rad(68), deg2rad(111)),
    Species(QE, ME,
            z->waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
            z->electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT)),
    Ground(15, 0.001),
    Transmitter(24e3)
)

const θs = [complex(r,i) for r = range(deg2rad(30), deg2rad(89), length=100) for i = range(deg2rad(-30), deg2rad(0), length=100)]

@testset "LongwaveModeSolver" begin
    include("test_EigenAngles.jl")
    include("test_Geophysics.jl")
    include("test_Waveguides.jl")

    include("test_magnetoionic.jl")
    include("test_wavefields.jl")
    include("test_modefinder.jl")
    include("test_modefinderX.jl")

    include("test_IO.jl")
end
