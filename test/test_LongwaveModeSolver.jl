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


function test_segmented(scenario)
    @unpack tx, rx, bfield, species, ground = scenario

    waveguide = LWMS.SegmentedWaveguide(LWMS.HomogeneousWaveguide)
    for i in eachindex(bfield)
        wvg = LWMS.HomogeneousWaveguide(bfield[i], species[i], ground[i])
        push!(waveguide, wvg)
    end

    E, amp, phase = bpm(waveguide, tx, rx)
end


function ttt(scenario)
    @unpack tx, rx, bfield, species, ground = scenario

    # b = bfield[1]
    # s = species[1]
    # g = ground[1]

    # jnk =  LWMS.HomogeneousWaveguide(b, s, g)

    waveguide = LWMS.SegmentedWaveguide(LWMS.HomogeneousWaveguide)
    for i in eachindex(bfield)
        wvg = LWMS.HomogeneousWaveguide(bfield[i], species[i], ground[i])
        push!(waveguide, wvg)
    end

    zs = range(LWMS.TOPHEIGHT, 0, length=513)
    nrsgmnt = length(waveguide)

    # Predetermine types
    Mtype = eltype(LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, waveguide[1]))
    wftype = promote_type(Mtype, Float64)
    ztype = typeof(zs)

    wavefields_vec = Vector{LWMS.Wavefields{wftype,ztype}}(undef, nrsgmnt)
    adjwavefields_vec = Vector{LWMS.Wavefields{wftype,ztype}}(undef, nrsgmnt)

    return nothing
end

    origcoords = LWMS.defaultcoordinates(tx.frequency)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = LWMS.GRPFParams(est_num_nodes, 1e-6, true)




end


    nsgmnt = 1
    wvg = waveguide[nsgmnt]

    interpolateM = true
    if interpolateM
        # TODO: check if functionwrapper is necessary
        Mfcn = LWMS.susceptibilityinterpolator(tx.frequency, wvg)
    else
        Mfcn = alt -> LWMS.susceptibility(alt, tx.frequency, wvg)
    end
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, wvg, Mfcn)

    modes = LWMS.findmodes(origcoords, grpfparams, modeequation)

    # adjoint wavefields are wavefields through adjoint waveguide, but for same modes
    # as wavefield
    @unpack bfield, species, ground = wvg
    adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
    adjwvg = LWMS.HomogeneousWaveguide(adjoint_bfield, species, ground)

    # TODO< just empty and resize the Wavefields
    wavefields = LWMS.Wavefields{wftype}(modes, zs)
    adjwavefields = LWMS.Wavefields{wftype}(modes, zs)

    LWMS.calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg)

    wavefields_vec[1] = wavefields
    adjwavefields_vec[1] = adjwavefields

    return nothing

    # return wavefields, adjwavefields, tx.frequency, wvg, adjwvg
end
