function verticalB()
    @unpack tx, rx, bfield, species, ground = verticalB_scenario()

    E, amp, phase = propagate(waveguide, tx, rx)

    return LMP.distance(rx, tx), amp, phase
end

function test_verticalB()
    lwpc_d, lwpc_a, lwpc_p = readlog("verticalB.log")
    lmp_d, lmp_a, lmp_p = verticalB()

    distmask = lwpc_d .> 400
    @test collect(lmp_d)/1e3 ≈ lwpc_d
    @test euclidean(lmp_a[distmask], lwpc_a[distmask]) < 3
    @test euclidean(rad2deg.(lmp_p[distmask]), lwpc_p[distmask]) < 15
end
