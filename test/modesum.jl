function test_modeterms(scenario)
    @unpack tx, rx, bfield, species, ground = scenario()
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    ea = TEST_ROOTS[scenario][1]
    modeequation = LMP.PhysicalModeEquation(ea, tx.frequency, waveguide)

    groundsampler = GroundSampler(rx.distance, rx.fieldcomponent)
    sampler = Sampler(rx.distance, 0.0, rx.fieldcomponent)

    tx1, rx1 = modeterms(modeequation, tx, sampler)
    tx2, rx2 = modeterms(modeequation, tx, groundsampler)

    @test tx1 ≈ tx2
    @test rx1 ≈ rx2
end
