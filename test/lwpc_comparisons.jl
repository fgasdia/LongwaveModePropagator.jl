function verticalB()
    @unpack bfield, species, ground = verticalB_scenario()

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    tx = Transmitter(VerticalDipole(), 24e3, 1e3)
    rx = GroundSampler(0:5e3:2000e3, Fields.Ex)

    E, amp, phase = propagate(waveguide, tx, rx)

    return LMP.distance(rx, tx), amp, phase
end

function test_verticalB()
    lwpc_d, lwpc_a, lwpc_p = readlog("verticalB.log")
    lmp_d, lmp_a, lmp_p = verticalB()

    @test lwpc_d ≈ collect(lmp_d)/1e3
end
