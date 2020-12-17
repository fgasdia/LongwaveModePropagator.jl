function test_modeconversion_segmented(scenario)
    @unpack distances, ea, tx, bfield, species, ground = scenario()
    params = LMPParams()

    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i],
                                        ground[i], distances[i]) for i in 1:2])


    heighttype = typeof(params.wavefieldheights)
    wavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, 2)
    adjwavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, 2)

    coordgrid = LMP.defaultcoordinates(tx.frequency)
    for i = 1:2
        wvg = waveguide[i]
        adjwvg = LMP.adjoint(wvg)
        modeequation = PhysicalModeEquation(tx.frequency, wvg)
        modes = findmodes(modeequation, coordgrid, params=params)

        wavefields = LMP.Wavefields(params.wavefieldheights, modes)
        adjwavefields = LMP.Wavefields(params.wavefieldheights, modes)
        LMP.calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg)

        wavefields_vec[i] = wavefields
        adjwavefields_vec[i] = adjwavefields
    end

    a = LMP.modeconversion(wavefields_vec[1], wavefields_vec[1], adjwavefields_vec[1])
    @test all(diag(a) .≈ complex(1))
    # a = modeconversion(wavefields_vec[1], wavefields_vec[2], adjwavefields_vec[2])
end
