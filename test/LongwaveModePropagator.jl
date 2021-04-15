function test_propagate(scenario)
    @unpack tx, rx, bfield, species, ground = scenario

    waveguide = HomogeneousWaveguide(bfield, species, ground)

    E, amp, phase = propagate(waveguide, tx, rx)
    @test eltype(E) == ComplexF64
    @test eltype(amp) == Float64
    @test eltype(phase) == Float64

    mesh = LMP.defaultmesh(tx.frequency)
    E2, amp2, phase2 = propagate(waveguide, tx, rx; mesh=mesh)
    @test E2 ≈ E    rtol=1e-3
    @test amp2 ≈ amp    rtol=1e-3
    @test phase2 ≈ phase    rtol=1e-3

    me = PhysicalModeEquation(tx.frequency, waveguide)
    modes = findmodes(me, mesh)
    E3, amp3, phase3 = propagate(waveguide, tx, rx; modes=modes)
    @test E3 ≈ E2    rtol=1e-3
    @test amp3 ≈ amp2    rtol=1e-3
    @test phase3 ≈ phase2    rtol=1e-3

    # mesh should be ignored
    E4, amp4, phase4 = propagate(waveguide, tx, rx; modes=modes, mesh=[1.0])
    @test E4 ≈ E3    rtol=1e-3
    @test amp4 ≈ amp3    rtol=1e-3
    @test phase4 ≈ phase3    rtol=1e-3

    # Are params being carried through?
    params = LMPParams(earthradius=6350e3)
    E5, amp5, phase5 = propagate(waveguide, tx, rx; params=params)
    @test !isapprox(E5, E, rtol=1e-3)
    @test !isapprox(amp5, amp, rtol=1e-3)
    @test !isapprox(phase5, phase, rtol=1e-3)

    # Mismatch between modes and modeterms
    E6, amp6, phase6 = propagate(waveguide, tx, rx; modes=modes, params=params)
    @test !isapprox(E6, E5, rtol=1e-3)
    @test !isapprox(amp6, amp5, rtol=1e-3)
    @test !isapprox(phase6, phase5, rtol=1e-3)
end

function test_propagate_segmented(scenario)
    @unpack tx, rx, bfield, species, ground, distances = scenario

    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i], ground[i],
                                                         distances[i]) for i in 1:2])

    E, amp, phase = propagate(waveguide, tx, rx)
    @test eltype(E) == ComplexF64
    @test eltype(amp) == Float64
    @test eltype(phase) == Float64

    mesh = LMP.defaultmesh(tx.frequency)
    E1, amp1, phase1 = propagate(waveguide, tx, rx)
    @test E1 ≈ E    rtol=1e-3
    @test amp1 ≈ amp    rtol=1e-3
    @test phase1 ≈ phase   rtol=1e-3

    # Are params being carried through?
    params = LMPParams(earthradius=6350e3)
    E3, amp3, phase3 = propagate(waveguide, tx, rx; params=params)
    @test !isapprox(E3, E, rtol=1e-3)
    @test !isapprox(amp3, amp, rtol=1e-3)
    @test !isapprox(phase3, phase, rtol=1e-3)
end

function test_mcranges_segmented(scenario)
    # Test that single receiver location has correct phase without continuous Range
    # This also confirms that results are correct if mode conversion occurs at a
    # segment_range that is not also an output_range

    @unpack tx, rx, bfield, species, ground, distances = scenario

    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i], ground[i],
                                                         distances[i]) for i in 1:2])

    Eref, ampref, phaseref = propagate(waveguide, tx, rx; unwrap=false)

    E1, a1, p1 = propagate(waveguide, tx, GroundSampler(600e3, Fields.Ez))
    E2, a2, p2 = propagate(waveguide, tx, GroundSampler(1400e3, Fields.Ez))
    E3, a3, p3 = propagate(waveguide, tx, GroundSampler(1800e3, Fields.Ez))

    m1 = findfirst(isequal(600e3), rx.distance)
    m2 = findfirst(isequal(1400e3), rx.distance)
    m3 = findfirst(isequal(1800e3), rx.distance)

    @test a1 ≈ ampref[m1] atol=0.01
    @test a2 ≈ ampref[m2] atol=0.01
    @test a3 ≈ ampref[m3] atol=0.01
    @test p1 ≈ phaseref[m1] atol=1e-3
    @test p2 ≈ phaseref[m2] atol=1e-3
    @test p3 ≈ phaseref[m3] atol=1e-3
end

@testset "LongwaveModePropagator.jl" begin
    @info "Testing LongwaveModePropagator"

    @info "  Running:"

    # Just to save time, running with only one scenario
    @info "    Homogeneous ionospheres..."
    for scn in (resonant_scenario, interp_scenario)
        test_propagate(scn)
    end

    @info "    Segmented ionospheres..."
    for scn in (segmented_scenario,)
        test_propagate_segmented(scn)
    end

    # propagate(file,...) is tested in IO.jl
end
