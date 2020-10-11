using Test
using Dates
using LinearAlgebra, Statistics
using StaticArrays
using Parameters
using OrdinaryDiffEq, DiffEqCallbacks
using RootsAndPoles, VoronoiDelaunay
using JSON3, StructTypes
using Interpolations, NLsolve, FiniteDiff

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const QE = LWMS.QE
const ME = LWMS.ME

# To profile in Juno
# @profiler (for i = 1:1000; LWMS.fcn(); end)

struct TestScenario{T1,T2,T3,T4,T5,T6}
    ea::T1
    bfield::T2
    species::T3
    ground::T4
    tx::T5
    rx::T6
end

const verticalB_scenario = TestScenario(
    EigenAngle(1.5 - 0.1im),
    BField(50e-6, π/2, 0),
    Species(QE, ME,
            z->waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
            z->electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT)),
    Ground(15, 0.001),
    Transmitter(24e3),
    GroundSampler(2000e3, Fields.Ez)
)

const resonant_scenario = TestScenario(
    EigenAngle(1.4161252139020892 - 0.016348911573820547im),  # resonant
    BField(50e-6, deg2rad(68), deg2rad(111)),
    Species(QE, ME,
            z->waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
            z->electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT)),
    Ground(15, 0.001),
    Transmitter(24e3),
    GroundSampler(2000e3, Fields.Ez)
)

const nonresonant_scenario = TestScenario(
    EigenAngle(1.5 - 0.1im),
    BField(50e-6, deg2rad(68), deg2rad(111)),
    Species(QE, ME,
            z->waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
            z->electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT)),
    Ground(15, 0.001),
    Transmitter(24e3),
    GroundSampler(2000e3, Fields.Ez)
)

const homogeneousiono_scenario = TestScenario(
    EigenAngle(1.4161252139020892 - 0.016348911573820547im),  # resonant
    BField(50e-6, deg2rad(68), deg2rad(111)),
    Species(QE, ME,
            z->z >= LWMS.CURVATURE_HEIGHT ? 2.65e6 : 0.0,
            z->z >= LWMS.CURVATURE_HEIGHT ? 1e8 : 0.0),
    Ground(15, 0.001),
    Transmitter(24e3),
    GroundSampler(2000e3, Fields.Ez)
)

const segmented_scenario = TestScenario(
    EigenAngle(1.5 - 0.1),
    [BField(50e-6, deg2rad(68), deg2rad(111)), BField(50e-6, deg2rad(68), deg2rad(111))],
    [Species(QE, ME,
             z->waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
             z->electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT)),
     Species(QE, ME,
             z->waitprofile(z, 77, 0.35, cutoff_low=LWMS.CURVATURE_HEIGHT),
             z->electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT))],
    [Ground(15, 0.001), Ground(15, 0.001)],
    Transmitter(24e3),
    GroundSampler(2000e3, Fields.Ez)
)

const θs = [complex(r,i) for r = range(deg2rad(30), deg2rad(89), length=100) for i = range(deg2rad(-30), deg2rad(0), length=100)]

@testset "LongwaveModeSolver" begin
    include("test_EigenAngles.jl")
    include("test_Geophysics.jl")
    include("test_Waveguides.jl")

    include("test_magnetoionic.jl")
    include("test_wavefields.jl")
    include("test_modefinder.jl")

    include("test_IO.jl")
end
