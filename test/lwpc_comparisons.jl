function verticalB()
    @unpack tx, rx, bfield, species, ground = verticalB_scenario()

    # params = LMPParams()
    # waveguide = HomogeneousWaveguide(bfield, species, ground)
    #
    # modeequation = LMP.PhysicalModeEquation(tx.frequency, waveguide)
    #
    # coordgrid = LMP.defaultcoordinates(tx.frequency)
    # modes = findmodes(modeequation, coordgrid, params=params)
    # 
    # return modes, tx, rx, waveguide

    E, amp, phase = propagate(waveguide, tx, rx)

    return LMP.distance(rx, tx), amp, phase
end

function temp(modes, tx, rx, waveguide)
    X = LMP.distance(rx, tx)
    params = LMPParams()

    E = zero(ComplexF64)

    txpower = LMP.power(tx)
    frequency = tx.frequency
    k = frequency.k

    i = 1
    ea = modes[i]
    modeequation = LMP.PhysicalModeEquation(ea, frequency, waveguide)
    txterm, rxterm = LMP.modeterms(modeequation, tx, rx, params=params)

    S₀ = LMP.referencetoground(ea.sinθ, params=params)
    expterm = -k*(S₀ - 1)
    txrxterm = txterm*rxterm

    E += txrxterm*cis(expterm*X[i])

    Q = 0.6822408*sqrt(frequency.f*txpower)  # factor from lw_sum_modes.for

    # TODO: Radiation resistance correction if zt > 0
    # See, e.g. Pappert Hitney 1989 TWIRE paper

    E *= Q/sqrt(abs(sin(X[i]/params.earthradius)))

    return txterm, rxterm, E
end

function test_verticalB()
    lwpc_d, lwpc_a, lwpc_p = readlog("verticalB.log")
    lmp_d, lmp_a, lmp_p = verticalB()

    @test lwpc_d ≈ collect(lmp_d)/1e3
end
