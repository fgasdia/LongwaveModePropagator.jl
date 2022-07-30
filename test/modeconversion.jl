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
    
    # Check that `a` is approximately identity. I'm not sure why the off-diagonal terms
    # aren't smaller
    di = diagind(a)
    for i in eachindex(a)
        if i in di
            @test a[i] ≈ 1
        else
            @test a[i] ≈ 0 atol=1e-2
        end
    end
    # a = modeconversion(wavefields_vec[1], wavefields_vec[2], adjwavefields_vec[2])
end

@testset "modeconversion.jl" begin
    @info "Testing modeconversion"

    for scn in (verticalB_scenario, multiplespecies_scenario)
        
        test_modeconversion_segmented(scn)
    end
end
