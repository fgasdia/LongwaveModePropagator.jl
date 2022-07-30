function test_modeconversion_segmented(scenario)
    @unpack tx, bfield, species, ground = scenario
    params = LMPParams()

    distances = (0.0, 500e3)
    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield, species, ground, distances[i]) for i in 1:2])

    heighttype = typeof(params.wavefieldheights)
    wavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, 2)
    adjwavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, 2)

    mesh = LMP.defaultmesh(tx.frequency)
    for i = 1:2
        wvg = waveguide[i]
        adjwvg = LMP.adjoint(wvg)
        modeequation = PhysicalModeEquation(tx.frequency, wvg)
        modes = findmodes(modeequation, mesh; params=params)

        wavefields = LMP.Wavefields(params.wavefieldheights, modes)
        adjwavefields = LMP.Wavefields(params.wavefieldheights, modes)
        LMP.calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg)

        wavefields_vec[i] = wavefields
        adjwavefields_vec[i] = adjwavefields
    end

    @time a = LMP.modeconversion(wavefields_vec[1], wavefields_vec[1], adjwavefields_vec[1])
    
    # Check that `a` is approximately identity.
    di = diagind(a)
    for i in eachindex(a)
        if i in di
            @test a[i] ≈ 1
        else
            @test a[i] ≈ 0 atol=1e-2
        end
    end

    # I'm not sure why the off-diagonal terms aren't smaller, but it doesn't appear to
    # result in any issues in the calculated electric field.
    # bfield = BField(50e-6, deg2rad(68), deg2rad(11))
    # tx = Transmitter(50e3)
    # rx = GroundSampler(0:1e3:2000e3, Fields.Ez)
    # ground = Ground(15, 0.001)
    # species = Species(QE, ME, z->waitprofile(z, 85, 0.5; cutoff_low=40e3), electroncollisionfrequency)
    # distances = (0.0, 1000e3, 1200e3, 1400e3)
    # waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield, species, ground, distances[i]) for i in 1:4])
    # E, a, p = propagate(waveguide, tx, rx)
    # plot(rx.distance/1e3, a, size=(900,600))
end

@testset "modeconversion.jl" begin
    @info "Testing modeconversion"

    for scn in (verticalB_scenario, multiplespecies_scenario)
        
        test_modeconversion_segmented(scn)
    end
end
