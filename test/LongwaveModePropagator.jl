function test_propagate(scenario)
    @unpack tx, rx, bfield, species, ground = scenario

    waveguide = HomogeneousWaveguide(bfield, species, ground)

    E, amp, phase = propagate(waveguide, tx, rx)
    @test eltype(E) == ComplexF64
    @test eltype(amp) == Float64
    @test eltype(phase) == Float64

    E1, amp1, phase1 = propagate(waveguide, tx, rx)
    @test meanabsdiff(E1, E) < 5
    @test meanabsdiff(amp1, amp) < 0.1
    @test meanabsdiff(phase1, phase) < 0.01

    mesh = LMP.defaultmesh(tx.frequency)
    E2, amp2, phase2 = propagate(waveguide, tx, rx; mesh=mesh)
    @test meanabsdiff(E2, E) < 5
    @test meanabsdiff(amp2, amp) < 0.1
    @test meanabsdiff(phase2, phase) < 0.01

    me = PhysicalModeEquation(tx.frequency, waveguide)
    modes = findmodes(me, mesh)
    E3, amp3, phase3 = propagate(waveguide, tx, rx; modes=modes)
    @test meanabsdiff(E3, E2) < 5
    @test meanabsdiff(amp3, amp2) < 0.1
    @test meanabsdiff(phase3, phase2) < 0.01

    # mesh should be ignored
    E4, amp4, phase4 = propagate(waveguide, tx, rx; modes=modes, mesh=[1.0])
    @test meanabsdiff(E4, E3) < 5
    @test meanabsdiff(amp4, amp3) < 0.1
    @test meanabsdiff(phase4, phase3) < 0.01

    # Multiple fields
    fullrx = GroundSampler(rx.distance, Fields.E)
    E7, amp7, phase7 = propagate(waveguide, tx, fullrx)
    @test meanabsdiff(E7[:,LMP.index(rx.fieldcomponent)], E) < 5
    @test meanabsdiff(amp7[:,LMP.index(rx.fieldcomponent)], amp) < 0.1
    @test meanabsdiff(phase7[:,LMP.index(rx.fieldcomponent)], phase) < 0.01

    distidx = findfirst(x->x==1000e3, rx.distance)
    fullrxpt = GroundSampler(1000e3, Fields.E)
    E8, amp8, phase8 = propagate(waveguide, tx, fullrxpt)
    @test meanabsdiff(E8[LMP.index(rx.fieldcomponent)], E[distidx]) < 5
    @test meanabsdiff(amp8[LMP.index(rx.fieldcomponent)], amp[distidx]) < 0.1
    @test meanabsdiff(mod2pi(phase8[LMP.index(rx.fieldcomponent)]), mod2pi(phase[distidx])) < 0.01

    rxpt = GroundSampler(1000e3, Fields.Ez)
    E9, amp9, phase9 = propagate(waveguide, tx, rxpt)
    @test meanabsdiff(E9, E[distidx]) < 5
    @test meanabsdiff(amp9, amp[distidx]) < 0.1
    @test meanabsdiff(mod2pi(phase9), mod2pi(phase[distidx])) < 0.01
end

function test_propagate_segmented(scenario)
    @unpack tx, rx, bfield, species, ground, distances = scenario

    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i], ground[i],
                                                         distances[i]) for i in 1:2])

    E, amp, phase = propagate(waveguide, tx, rx)
    @test eltype(E) == ComplexF64
    @test eltype(amp) == Float64
    @test eltype(phase) == Float64

    E1, amp1, phase1 = propagate(waveguide, tx, rx)
    @test meanabsdiff(E1, E) < 1
    @test meanabsdiff(amp1, amp) < 1
    @test meanabsdiff(phase1, phase) < 0.1
end

function test_mcranges_segmented(scenario)
    # Test that single receiver location has correct phase without continuous Range
    # This also confirms that results are correct if mode conversion occurs at a
    # segment_range that is not also an output_range

    @unpack tx, rx, bfield, species, ground, distances = scenario

    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i], ground[i],
                                                         distances[i]) for i in 1:2])

    _, ampref, phaseref = propagate(waveguide, tx, rx; unwrap=false)

    _, a1, p1 = propagate(waveguide, tx, GroundSampler(600e3, Fields.Ez))
    _, a2, p2 = propagate(waveguide, tx, GroundSampler(1400e3, Fields.Ez))
    _, a3, p3 = propagate(waveguide, tx, GroundSampler(1800e3, Fields.Ez))

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
