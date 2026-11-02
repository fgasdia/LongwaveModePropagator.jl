function test_propagate(scenario)
    @unpack tx, rx, bfield, species, ground = scenario

    waveguide = HomogeneousWaveguide(bfield, species, ground)

    E, amp, phase = propagate(waveguide, tx, rx)
    @test eltype(E) == ComplexF64
    @test eltype(amp) == Float64
    @test eltype(phase) == Float64

    E1, amp1, phase1 = propagate(waveguide, tx, rx)
    @test E1 ≈ E    atol=500 nans=true
    @test amp1 ≈ amp    atol=0.5 nans=true
    @test phase1 ≈ phase    atol=0.05 nans=true

    mesh = LMP.defaultmesh(tx.frequency)
    E2, amp2, phase2 = propagate(waveguide, tx, rx; mesh=mesh)
    @test E2 ≈ E    atol=500 nans=true
    @test amp2 ≈ amp    atol=0.5 nans=true
    @test phase2 ≈ phase    atol=0.05 nans=true

    me = PhysicalModeEquation(tx.frequency, waveguide)
    modes = findmodes(me, mesh)
    E3, amp3, phase3 = propagate(waveguide, tx, rx; modes=modes)
    @test E3 ≈ E2    atol=500 nans=true
    @test amp3 ≈ amp2    atol=0.5 nans=true
    @test phase3 ≈ phase2    atol=0.05 nans=true

    # mesh should be ignored
    E4, amp4, phase4 = propagate(waveguide, tx, rx; modes=modes, mesh=[1.0])
    @test E4 ≈ E3    atol=500 nans=true
    @test amp4 ≈ amp3    atol=0.5 nans=true
    @test phase4 ≈ phase3    atol=0.05 nans=true

    # Multiple fields
    fullrx = GroundSampler(rx.distance, Fields.E)
    E7, amp7, phase7 = propagate(waveguide, tx, fullrx)
    @test E7[:,LMP.index(rx.fieldcomponent)] ≈ E    atol=500 nans=true
    @test amp7[:,LMP.index(rx.fieldcomponent)] ≈ amp    atol=0.5 nans=true
    @test phase7[:,LMP.index(rx.fieldcomponent)] ≈ phase    atol=0.1 nans=true

    distidx = findfirst(x->x==1000e3, rx.distance)
    fullrxpt = GroundSampler(1000e3, Fields.E)
    E8, amp8, phase8 = propagate(waveguide, tx, fullrxpt)
    @test E8[LMP.index(rx.fieldcomponent)] ≈ E[distidx]    atol=500 nans=true
    @test amp8[LMP.index(rx.fieldcomponent)] ≈ amp[distidx]    atol=0.5 nans=true
    @test mod2pi(phase8[LMP.index(rx.fieldcomponent)]) ≈ phase[distidx]    atol=0.05 nans=true

    rxpt = GroundSampler(1000e3, Fields.Ez)
    E9, amp9, phase9 = propagate(waveguide, tx, rxpt)
    @test E9 ≈ E[distidx]    atol=500 nans=true
    @test amp9 ≈ amp[distidx]    atol=0.5 nans=true
    @test mod2pi(phase9) ≈ phase[distidx]    atol=0.05 nans=true
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
    @test E1 ≈ E    atol=500 nans=true
    @test amp1 ≈ amp    atol=0.5 nans=true
    @test phase1 ≈ phase    atol=0.05 nans=true
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
