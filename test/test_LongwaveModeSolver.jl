function bpm_test(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)

    E, phase, amp = bpm(waveguide, tx, rx)

    # Read LWPC results
    lwpcfile = "homogeneous_day_lwpc.json"  # TODO adapt to scenario
    if isfile(lwpcfile)
        lwpcres = open(lwpcfile,"r") do f
            v = JSON3.read(f, LWMS.BasicOutput)
            return v
        end
    else
        @info "$lwpcfile not found"
        return false
    end

    # origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    # est_num_nodes = ceil(Int, length(origcoords)*1.5)
    # grpfparams = LWMS.GRPFParams(est_num_nodes, 1e-6, true)
    #
    # Mfcn = LWMS.susceptibilityinterpolator(tx.frequency, waveguide)
    # modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
    #
    # modes = LWMS.findmodes(origcoords, grpfparams, modeequation)
    #
    # ea = modes[1]
    #
    # Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, waveguide)
    # modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
    #
    # dFdθ, R, Rg = LWMS.solvemodalequation(ea, modeequation, LWMS.Dθ())
    # # efconstants = LWMS.excitationfactorconstants(ea, R, Rg, tx.frequency, waveguide.ground)
    # #
    # # Sγ, Cγ = sincos(π/2 - inclination(tx))  # γ is measured from vertical
    # # Sϕ, Cϕ = sincos(azimuth(tx))  # ϕ is measured from `x`
    # #
    # # emitter_orientation = (t1=Cγ, t2=Sγ*Sϕ, t3=Sγ*Cϕ, zt=zt)
    # # sampler_orientation = (rxcomponent=rxcomponent, zr=zr)
    #
    # return ea, R, Rg, tx.frequency, waveguide.ground

    # return modes, waveguide, tx, rx

    # @bp E = Efield(modes, waveguide, tx, rx)
end

# excitationfactorconstants(e, R, Rg, f, g)
