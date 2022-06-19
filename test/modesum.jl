function test_modeterms(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    ea = TEST_MODES[scenario][1]
    modeequation = PhysicalModeEquation(ea, tx.frequency, waveguide)

    groundsampler = GroundSampler(rx.distance, rx.fieldcomponent)
    sampler = Sampler(rx.distance, rx.fieldcomponent, 0.0)

    tx1, rx1 = LMP.modeterms(modeequation, tx, sampler)
    tx2, rx2 = LMP.modeterms(modeequation, tx, groundsampler)  # specialized

    @test length(tx1) == 1
    @test length(rx1) == LMP.NUMRXTERMS
    
    @test tx1 ≈ tx2
    @test rx1 ≈ rx2

    # frequency mismatch with modeequation
    txwrong = Transmitter(15e3)
    @test_throws ArgumentError LMP.modeterms(modeequation, txwrong, sampler)
end

function test_Efield()
    @unpack bfield, species, ground, tx, rx = verticalB_scenario
    homogeneouswaveguide = HomogeneousWaveguide(bfield, species, ground)
    homogeneousmodes = TEST_MODES[verticalB_scenario]

    E1 = LMP.Efield(homogeneousmodes, homogeneouswaveguide, tx, rx)

    # Segmented waveguide where each segment is identical to the homogeneous waveguide
    segmentedwaveguide = SegmentedWaveguide([
        HomogeneousWaveguide(bfield, species, ground),
        HomogeneousWaveguide(bfield, species, ground),
        HomogeneousWaveguide(bfield, species, ground)
    
    ])
    J = length(segmentedwaveguide)

    params = LMPParams()
    # for testing impact on difference between wavefields, below
    # params = LMPParams(wavefieldintegrationparams=IntegrationParams(solver=Tsit5(), tolerance=1e-8))
    # params = LMPParams(grpfparams=GRPFParams(100000, 1e-6, true),
        # integrationparams=IntegrationParams(solver=Vern7(), tolerance=1e-5))
    tolerance = params.grpfparams.tolerance
    mesh = LMP.defaultmesh(tx.frequency)

    heighttype = typeof(params.wavefieldheights)
    wavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, J)
    adjwavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, J)

    # Calculate wavefields and adjoint wavefields for each segment of waveguide
    for j in 1:J
        wvg = segmentedwaveguide[j]

        modeequation = PhysicalModeEquation(tx.frequency, wvg)

        modes = findmodes(modeequation, mesh; params=params)
        adjwvg = LMP.adjoint(wvg)

        wavefields = LMP.Wavefields(params.wavefieldheights, modes)
        adjwavefields = LMP.Wavefields(params.wavefieldheights, modes)

        LMP.calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg;
                              params=params)

        wavefields_vec[j] = wavefields
        adjwavefields_vec[j] = adjwavefields
    end

    # TEMP: Move to test/wavefields?
    @test isapprox(getfield.(wavefields_vec[1].eas, :θ),
        getfield.(wavefields_vec[2].eas, :θ), atol=tolerance)
    @test isapprox(getfield.(wavefields_vec[1].eas, :θ),
        getfield.(wavefields_vec[3].eas, :θ), atol=tolerance)

    @test size(wavefields_vec[1].v) == size(wavefields_vec[2].v) == size(wavefields_vec[3].v)
    for i in eachindex(wavefields_vec[1].v)
        # if maxabsdiff(wavefields_vec[1].v[i], wavefields_vec[2].v[i]) > 1e-3
            # println(CartesianIndices(size(wavefields_vec[1].v))[i])
        # end
        @test maxabsdiff(wavefields_vec[1].v[i], wavefields_vec[2].v[i]) < 0.015
    end

    # examine the absolute difference in wavefields for "bad" modes
    # The largest differences appear where the magnitude of the wavefields is largest.
    # There is no odd structure in the difference plots, so I am not particularly concerned.
    # Decreasing the GRPF tolerance by one order of magnitude improves the wavefield match,
    # but at a significant hit to overall runtime (50%). Decreasing the wavefieldintegration
    # tolerance has a minor improvement and only costs a few tens of milliseconds
    # using Plots
    # plot()
    # for i = 1:6
    #     plot!(abs.(getindex.(wavefields_vec[1].v[:,4],i) .- getindex.(wavefields_vec[2].v[:,4],i)),
    #         params.wavefieldheights/1e3, label=i)
    # end
    # plot!()

    E2 = LMP.Efield(segmentedwaveguide, wavefields_vec, adjwavefields_vec, tx, rx; params=params)


    # E2 is better than 1.5% different from E1 - this is just at the measurement noise
    @test all(abs.((E2 .- E1) ./ E1) .< 0.015)
    @test maxabsdiff(LMP.dBamplitude.(E2), LMP.dBamplitude.(E1)) < 0.11
    @test rad2deg.(maxabsdiff(angle.(E2), angle.(E1))) < 1
end


@testset "modesum.jl" begin
    @info "Testing modesum"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        multiplespecies_scenario)
        
        test_modeterms(scn)
    end
    test_Efield()
end
