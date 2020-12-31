using Test, Dates, Random
using LinearAlgebra, Statistics
using StaticArrays
using Parameters
using OrdinaryDiffEq, DiffEqCallbacks
using JSON3, StructTypes, CSV
using Interpolations, NLsolve, FiniteDiff

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, distance, power
const LMP = LongwaveModePropagator

const TEST_RNG = MersenneTwister(1234)

const LWPC_PATH = "LWPC"

include("utils.jl")

#==
Scenarios
==#
const TEST_MODES = Dict{Any,Vector{EigenAngle}}()

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
    rx=GroundSampler(0:5e3:2000e3, Fields.Ez)
)

const resonant_scenario = (
    ea=EigenAngle(1.416127852502346 - 0.016482589477369265im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=GroundSampler(0:5e3:2000e3, Fields.Ez)
)

const resonant_elevatedrx_scenario = (
    ea=EigenAngle(1.416127852502346 - 0.016482589477369265im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=Sampler(0:5e3:2000e3, Fields.Ez, 8e3)
)

const resonant_horizontal_scenario = (
    ea=EigenAngle(1.416127852502346 - 0.016482589477369265im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=GroundSampler(0:5e3:2000e3, Fields.Ey)
)

const nonresonant_scenario = (
    ea=EigenAngle(1.5 - 0.1im),
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                    electroncollisionfrequency),
    ground=Ground(5, 0.00005),
    tx=Transmitter(24e3),
    rx=GroundSampler(0:5e3:2000e3, Fields.Ez)
)

const homogeneousiono_scenario = (
    ea=EigenAngle(1.4161252139020892 - 0.016348911573820547im),  # resonant
    bfield=BField(50e-6, deg2rad(68), deg2rad(111)),
    species=Species(QE, ME,
                    z->2.65e6,
                    z->1e8),
    ground=Ground(15, 0.001),
    tx=Transmitter(24e3),
    rx=GroundSampler(0:5e3:2000e3, Fields.Ez)
)

const segmented_scenario = (
    ea=EigenAngle(1.5 - 0.1),
    bfield=[BField(50e-6, deg2rad(68), deg2rad(111)), BField(50e-6, deg2rad(68), deg2rad(111))],
    species=[Species(QE, ME,
                     z->waitprofile(z, 75, 0.32, cutoff_low=40e3),
                     electroncollisionfrequency),
             Species(QE, ME,
                     z->waitprofile(z, 80, 0.45, cutoff_low=40e3),
                     electroncollisionfrequency)],
    ground=[Ground(15, 0.001), Ground(15, 0.001)],
    tx=Transmitter(24e3),
    rx=GroundSampler(0:5e3:2000e3, Fields.Ez),
    distances=[0.0, 1000e3],
)

const θs = [complex(r,i) for r = range(deg2rad(60), deg2rad(89), length=30) for i =
            range(deg2rad(-10), deg2rad(0), length=11)]

maxabsdiff(a, b) = maximum(abs.(a - b))
meanabsdiff(a, b) = mean(abs.(a - b))

function findroots(scenario)
    @unpack bfield, species, ground, tx = scenario
    w = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(tx.frequency, w)
    origcoords = LMP.defaultmesh(tx.frequency)
    return findmodes(me, origcoords)
end


@testset "LongwaveModePropagator" begin
    include("EigenAngles.jl")
    include("Geophysics.jl")
    include("Waveguides.jl")

    include("magnetoionic.jl")
    include("TMatrix.jl")
    include("bookerquartic.jl")

    include("modefinder.jl")
    include("wavefields.jl")
    include("modeconversion.jl")
    include("modesum.jl")

    include("LongwaveModePropagator.jl")
    include("lwpc_comparisons.jl")

    include("IO.jl")
end
