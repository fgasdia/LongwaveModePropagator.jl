using Test, Dates, Random
using LinearAlgebra, Statistics
using StaticArrays
using Parameters
using OrdinaryDiffEq, DiffEqCallbacks
using JSON3, StructTypes, CSV
using Interpolations, NLsolve, FiniteDiff

using LongwaveModePropagator
const LMP = LongwaveModePropagator

const QE = LMP.QE
const ME = LMP.ME

const TEST_RNG = MersenneTwister(1234)


const verticalB_scenario = (
    ea=EigenAngle(1.5 - 0.1im),
    bfield=BField(50e-6, π/2, 0),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=GroundSampler(0:5e3:2000e3, Fields.Ez)
)

const isotropicB_resonant_scenario = (
    ea=EigenAngle(1.453098822238508 - 0.042008075239068944im),  # resonant
    bfield=BField(50e-6, 0, π/2),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=GroundSampler(2000e3, Fields.Ez)
)

const resonant_scenario = (
    ea=EigenAngle(1.416127852502346 - 0.016482589477369265im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=GroundSampler(2000e3, Fields.Ez)
)

const nonresonant_scenario = (
    ea=EigenAngle(1.5 - 0.1im),
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(10, 0.0003),
    tx=Transmitter(24e3),
    rx=GroundSampler(2000e3, Fields.Ez)
)

const homogeneousiono_scenario = (
    ea=EigenAngle(1.4161252139020892 - 0.016348911573820547im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->2.65e6,
                    z->1e8),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=GroundSampler(2000e3, Fields.Ez)
)

const segmented_scenario = (
    ea=EigenAngle(1.5 - 0.1),
    bfield=[BField(50e-6, deg2rad(68), deg2rad(111)), BField(50e-6, deg2rad(68), deg2rad(111))],
    species=[Species(QE, ME,
                     z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                     electroncollisionfrequency),
             Species(QE, ME,
                     z->waitprofile(z, 77, 0.35, cutoff_low=40e3),
                     electroncollisionfrequency)],
    ground=[Ground(15, 0.001), Ground(15, 0.001)],
    tx=Transmitter(24e3),
    rx=GroundSampler(2000e3, Fields.Ez),
    distances=[0.0, 1000e3],
)

const θs = [complex(r,i) for r = range(deg2rad(60), deg2rad(89), length=30) for i = range(deg2rad(-15), deg2rad(0), length=16)]
err_func(a,b) = maximum(abs.(a-b))

function findroots(scenario)
    @unpack bfield, species, ground, tx = scenario
    w = LMP.HomogeneousWaveguide(bfield, species, ground)
    me = LMP.PhysicalModeEquation(tx.frequency, w)
    origcoords = LMP.defaultcoordinates(tx.frequency)
    return LMP.findmodes(me, origcoords)
end

@testset "LongwaveModePropagator" begin
    include("LWPC_utils.jl")

    include("EigenAngles.jl")
    include("Geophysics.jl")
    include("Waveguides.jl")

    include("magnetoionic.jl")
    include("TMatrix.jl")
    include("bookerquartic.jl")

    @test_nowarn const TEST_ROOTS = Dict{Any,Vector{EigenAngle}}((scenario, findroots(scenario))
        for scenario in (verticalB_scenario, resonant_scenario, nonresonant_scenario))

    include("modefinder.jl")
    include("wavefields.jl")
    include("modeconversion.jl")

    include("test_IO.jl")
end
