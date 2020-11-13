function scenario()
    ea = EigenAngle(deg2rad(complex(85.0, -1.0)))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    # `tx` isn't just bits
    ground = Ground(15, 0.001)
    bfield = BField(50e-6, π/2, 0)

    electrons = Species(QE, ME,
                        z -> waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
                        z -> electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT))

    return ea, tx, ground, bfield, electrons
end

function dXdzW(X, dzparams, z)
    # dX/dz using the intermediate matrix W
    @unpack ea, frequency, Mfcn = dzparams

    k = frequency.k
    C = ea.cosθ

    M = Mfcn(z)
    T = LWMS.tmatrix(ea, M)
    W11, W21, W12, W22 = LWMS.wmatrix(ea, T)

    a = (W21 - W22 + W11 - W12)/C
    b = W22 + W12
    c = -W11 + W12
    d = -C*W12

    return -1im/2*k*(a + b*X + X*c + X*d*X)
end

function test_dX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, electrons)

    Mtop = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, electrons)
    Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
    Xtop = LWMS.R2X(ea, Rtop)

    dzparams = LWMS.DZParams(ea, tx.frequency, Mfcn)
    X = LWMS.dXdz(Xtop, dzparams, LWMS.TOPHEIGHT-500)
    Xref = dXdzW(Xtop, dzparams, LWMS.TOPHEIGHT-500)

    return X ≈ Xref
end

function test_R2XX2R()
    ea = EigenAngle(deg2rad(85.0-1.0im))
    R = rand(SMatrix{2,2,ComplexF64})
    X = LWMS.R2X(ea, R)
    Rtest = LWMS.X2R(ea, X)

    return Rtest ≈ R
end

function test_modalequationX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, electrons)
    modeequation = LWMS.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)

    # known solution
    ea = EigenAngle(1.4152764714690873 - 0.01755942938376613im)

    f = LWMS.solvemodalequation(ea, modeequation)

    return LWMS.isroot(f)
end

function test_modefinderX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, electrons)
    modeequation = LWMS.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)

    # Δr from 0.5->0.25 => time from 3.8->5.3 sec
    # tolerance from 1e-8->1e-7 => time from 5.3->4.6 sec
    origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = LWMS.GRPFParams(est_num_nodes, 1e-5, true)

    modes = LWMS.findmodes(origcoords, grpfparams, modeequation)

    physicalmodeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
    for m in modes
        f = LWMS.solvemodalequation(m, physicalmodeequation)
        LWMS.isroot(f) || return false
    end
    return true
end


function integratedreflectionX_deriv(scenario)
    @unpack tx, ground, bfield, species = scenario
    freq = tx.frequency

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.ModifiedModeEquation(freq, waveguide)

    Cs = [cos(θ) for θ in θs]

    @info "    Integrating dX/dC for many eigenangles. This may take a minute."

    Xref(θ,m) = (m = LWMS.setea(EigenAngle(θ), m); LWMS.integratedreflectionX(m))
    XdX(θ,m) = (m = LWMS.setea(EigenAngle(θ), m); LWMS.integratedreflectionX(m, LWMS.DC()))

    Xs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(θs))
    dXs = similar(Xs)
    Xrefs = similar(Xs)
    dXrefs = similar(Xs)
    Threads.@threads for i in 1:length(Cs)  # NOTE: this also effectively checks for thread safety
        v = XdX(θs[i], modeequation)
        Xs[i] = v[SVector(1,2),:]
        dXs[i] = v[SVector(3,4),:]
        Xrefs[i] = Xref(θs[i], modeequation)
        dXrefs[i] = FiniteDiff.finite_difference_derivative(z->Xref(z, modeequation), θs[i],
                                                            Val{:central})
    end

    for i = 1:4
        R = [v[i] for v in Rs]
        dR = [v[i] for v in dRs]
        Rr = [v[i] for v in Rrefs]
        dRr = [v[i] for v in dRrefs]

        # Rref.(θs) !≈ R.(θs) because of difference in integration tolerance.
        # For very large Rs, the difference can be significant, therefore we use rtol
        isapprox(R, Rr, rtol=1e-2) || return false
        isapprox(dR, dRr, rtol=1e-2) || return false
    end

    return true
end








function test_modefinderX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = GRPFParams(est_num_nodes, 1e-8, true)

    @unpack bfield, species = waveguide
    Mfcn = LWMS.susceptibilityinterpolator(tx.frequency, bfield, electrons)

    modeequation = LWMS.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)
    f(θ) = solvemodalequation(EigenAngle(θ), modeequation)
    roots, poles = grpf(f, origcoords, grpfparams)

    # Ensure roots are valid solutions to the modal equation
    filterroots!(roots, frequency, waveguide)

    # Remove any redundant modes
    # if tolerance is 1e-8, this rounds to 7 decimal places
    sort!(roots, by=reim)
    ndigits = round(Int, abs(log10(grpfparams.tolerance)+1), RoundDown)
    unique!(z->round(z, digits=ndigits), roots)
end

function test_integratedX(scenario)
    @unpack tx, ground, bfield, species = scenario

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.ModifiedModeEquation(EigenAngle(deg2rad(85-0.5im)), tx.frequency, waveguide)
    sol = LWMS.integratedreflectionX(modeequation)

    return sol
end

function test_freespace(scenario)
    @unpack tx, ground, bfield, species = scenario

    freq = tx.frequency
    ea = EigenAngle(deg2rad(85-0.5im))

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.ModifiedModeEquation(EigenAngle(deg2rad(85-0.5im)), tx.frequency, waveguide)
    Xb = LWMS.integratedreflectionX(modeequation)

    # TODO: Do a test over different ea to see how much the altitude varies
    # TODO: modefinder at different reference altitudes to determine sensitivity

    alts = 30e3:90e3
    Xs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(alts))
    for i in eachindex(Xs)
        X = LWMS.freespaceintegration(alts[i], Xb, ea, freq)
        Xs[i] = X
    end

    return Xs
end

aX1 = [abs2(x[1]) for x in Xs]
aX2 = [abs2(x[2]) for x in Xs]
aX3 = [abs2(x[3]) for x in Xs]
aX4 = [abs2(x[4]) for x in Xs]

function tmp(scenario)
    @unpack tx, ground, bfield, species = scenario
    @unpack tolerance, solver, force_dtmin = LWMS.DEFAULT_INTEGRATIONPARAMS

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.ModifiedModeEquation(EigenAngle(deg2rad(85-0.5im)), tx.frequency, waveguide)

    centerea = EigenAngle(deg2rad(70.0-4im))
    modeequation = LWMS.setea(centerea, modeequation)

    Xb = LWMS.integratedreflectionX(modeequation)

    mindXsum = Inf
    mindXaltitude = Inf

    altitudes = 40e3:1e3:92e3

    frequency = modeequation.frequency

    θs = SVector(deg2rad(60.0-0.1im), deg2rad(60.0-10im),
                 deg2rad(89.9-0.1im), deg2rad(89.9-10im))

    for i in eachindex(altitudes)
        f = θ -> LWMS.freespaceintegration(altitudes[i], Xb, EigenAngle(θ), frequency)
        dXs = FiniteDiff.finite_difference_derivative(f, θs, Val{:central}, typeof(Xb))

        dXsum = 0.0
        for t in 1:4
            dXsum += (abs2(dXs[t][1,1]) + abs2(dXs[t][2,1]) + abs2(dXs[t][1,2]) +
                      abs2(dXs[t][2,2]))
        end

        if dXsum < mindXsum
            mindXsum = dXsum
            mindXaltitude = altitudes[i]
        end
    end
    return mindXsum, mindXaltitude
end

function reflectionheight(scenario)
    @unpack tx, ground, bfield, species = scenario

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.ModifiedModeEquation(EigenAngle(deg2rad(85-0.5im)), tx.frequency, waveguide)
    #
    # minXsum, minXaltitude = LWMS.reflectionheight(modeequation)

    LWMS.reflectionheight(modeequation)
end

function tdx()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    dzparams = LWMS.DZParams(ea, tx.frequency, Mfcn)

    Xtop = rand(SMatrix{2,2,ComplexF64})
    # Xtop = @SMatrix [2.0+1.0im 2.1+1.1im
                     # 19+0.2im 13.0+0.5im]
    X = LWMS.dXdz(Xtop, dzparams, 80e3)
    # R = LWMS.dRdz(Xtop, dzparams, 80e3)
end

function tsx()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    # dzparams = LWMS.DZParams(ea, tx.frequency, Mfcn)

    modeequation = LWMS.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)
    LWMS.solvemodalequation(ea, modeequation)
end


@testset "modefinderX.jl" begin
    @info "Testing modefinderX"

    @test test_dX()
    @test test_R2XX2R()
    @test test_modalequationX()
end
