using Test, Dates, Random
using LinearAlgebra, Statistics
using StaticArrays
using Parameters
using OrdinaryDiffEq, DiffEqCallbacks
using RootsAndPoles, VoronoiDelaunay
using JSON3, StructTypes
using Interpolations, NLsolve, FiniteDiff

using LongwaveModePropagator
const LMP = LongwaveModePropagator

const QE = LMP.QE
const ME = LMP.ME

const TEST_RNG = MersenneTwister(1234)


const verticalB_scenario = @with_kw (
    ea=EigenAngle(1.5 - 0.1im),
    bfield=BField(50e-6, π/2, 0),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    # Transmitter(24e3),
    # GroundSampler(2000e3, Fields.Ez)
)

const isotropicB_resonant_scenario = @with_kw (
    ea=EigenAngle(1.453098822238508 - 0.042008075239068944im),  # resonant
    bfield=BField(50e-6, 0, π/2),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    # Transmitter(24e3),
    # GroundSampler(2000e3, Fields.Ez)
)

const resonant_scenario = @with_kw (
    ea=EigenAngle(1.416127852502346 - 0.016482589477369265im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    # Transmitter(24e3),
    # GroundSampler(2000e3, Fields.Ez)
)

const nonresonant_scenario = @with_kw (
    ea=EigenAngle(1.5 - 0.1im),
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    # Transmitter(24e3),
    # GroundSampler(2000e3, Fields.Ez)
)

const homogeneousiono_scenario = @with_kw (
    ea=EigenAngle(1.4161252139020892 - 0.016348911573820547im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->2.65e6,
                    z->1e8),
    ground=Ground(15, 0.001),
    # Transmitter(24e3),
    # GroundSampler(2000e3, Fields.Ez)
)

const segmented_scenario = @with_kw (
    ea=EigenAngle(1.5 - 0.1),
    bfield=[BField(50e-6, deg2rad(68), deg2rad(111)), BField(50e-6, deg2rad(68), deg2rad(111))],
    species=[Species(QE, ME,
                     z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                     electroncollisionfrequency),
             Species(QE, ME,
                     z->waitprofile(z, 77, 0.35, cutoff_low=40e3),
                     electroncollisionfrequency)],
    ground=[Ground(15, 0.001), Ground(15, 0.001)],
    # Transmitter(24e3),
    # GroundSampler(2000e3, Fields.Ez)
)

const θs = [complex(r,i) for r = range(deg2rad(60), deg2rad(89), length=30) for i = range(deg2rad(-15), deg2rad(0), length=16)]
err_func(a,b) = maximum(abs.(a-b))


@testset "LongwaveModePropagator" begin
    include("test_EigenAngles.jl")
    include("test_Geophysics.jl")
    include("test_Waveguides.jl")

    include("test_magnetoionic.jl")
    include("test_modefinder.jl")
    include("test_wavefields.jl")

    include("test_IO.jl")
end
