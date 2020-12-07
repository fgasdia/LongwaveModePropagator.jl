function scenario()
    ea = EigenAngle(deg2rad(complex(85.0, -1.0)))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    # `tx` isn't just bits
    ground = Ground(15, 0.001)
    bfield = BField(50e-6, π/2, 0)

    electrons = Species(QE, ME,
                        z -> waitprofile(z, 75, 0.32, cutoff_low=LMP.CURVATURE_HEIGHT),
                        z -> electroncollisionfrequency(z, cutoff_low=LMP.CURVATURE_HEIGHT))

    return ea, tx, ground, bfield, electrons
end

function dXdzW(X, dzparams, z)
    # dX/dz using the intermediate matrix W
    @unpack ea, frequency, Mfcn = dzparams

    k = frequency.k
    C = ea.cosθ

    M = Mfcn(z)
    T = LMP.tmatrix(ea, M)
    W11, W21, W12, W22 = LMP.wmatrix(ea, T)

    a = (W21 - W22 + W11 - W12)/C
    b = W22 + W12
    c = -W11 + W12
    d = -C*W12

    return -1im/2*k*(a + b*X + X*c + X*d*X)
end

function test_dX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(alt) = LMP.susceptibility(alt, tx.frequency, bfield, electrons)

    Mtop = LMP.susceptibility(LMP.TOPHEIGHT, tx.frequency, bfield, electrons)
    Rtop = LMP.sharpboundaryreflection(ea, Mtop)
    Xtop = LMP.R2X(ea, Rtop)

    dzparams = LMP.DZParams(ea, tx.frequency, Mfcn)
    X = LMP.dXdz(Xtop, dzparams, LMP.TOPHEIGHT-500)
    Xref = dXdzW(Xtop, dzparams, LMP.TOPHEIGHT-500)

    return X ≈ Xref
end

function test_R2XX2R()
    ea = EigenAngle(deg2rad(85.0-1.0im))
    R = rand(SMatrix{2,2,ComplexF64})
    X = LMP.R2X(ea, R)
    Rtest = LMP.X2R(ea, X)

    return Rtest ≈ R
end

function test_modalequationX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, electrons, ground)

    Mfcn(alt) = LMP.susceptibility(alt, tx.frequency, bfield, electrons)
    modeequation = LMP.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)

    # known solution
    ea = EigenAngle(1.4152764714690873 - 0.01755942938376613im)

    f = LMP.solvemodalequation(ea, modeequation)

    return LMP.isroot(f)
end

function test_modefinderX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, electrons, ground)

    Mfcn(alt) = LMP.susceptibility(alt, tx.frequency, bfield, electrons)
    modeequation = LMP.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)

    # Δr from 0.5->0.25 => time from 3.8->5.3 sec
    # tolerance from 1e-8->1e-7 => time from 5.3->4.6 sec
    origcoords = LMP.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = LMP.GRPFParams(est_num_nodes, 1e-5, true)

    modes = LMP.findmodes(origcoords, grpfparams, modeequation)

    physicalmodeequation = LMP.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
    for m in modes
        f = LMP.solvemodalequation(m, physicalmodeequation)
        LMP.isroot(f) || return false
    end
    return true
end










function test_modefinderX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LMP.susceptibility(z, tx.frequency, bfield, electrons)

    origcoords = LMP.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = GRPFParams(est_num_nodes, 1e-8, true)

    @unpack bfield, species = waveguide
    Mfcn = LMP.susceptibilityinterpolator(tx.frequency, bfield, electrons)

    modeequation = LMP.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)
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

function test_integratedX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LMP.susceptibility(z, tx.frequency, bfield, electrons)

    origcoords = LMP.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = GRPFParams(est_num_nodes, 1e-8, true)

    @unpack bfield, species = waveguide
    # Mfcn = LMP.susceptibilityinterpolator(tx.frequency, bfield, species)
    susceptibilitymfcn = alt->LMP.susceptibility(alt, tx.frequency, bfield, species)

    LMP.integratedreflectionX(ea, tx.frequency, waveguide, susceptibilitymfcn)
end

function tdx()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LMP.susceptibility(z, tx.frequency, bfield, electrons)

    dzparams = LMP.DZParams(ea, tx.frequency, Mfcn)

    Xtop = rand(SMatrix{2,2,ComplexF64})
    # Xtop = @SMatrix [2.0+1.0im 2.1+1.1im
                     # 19+0.2im 13.0+0.5im]
    X = LMP.dXdz(Xtop, dzparams, 80e3)
    # R = LMP.dRdz(Xtop, dzparams, 80e3)
end

function tsx()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LMP.susceptibility(z, tx.frequency, bfield, electrons)

    # dzparams = LMP.DZParams(ea, tx.frequency, Mfcn)

    modeequation = LMP.ModifiedModeEquation(tx.frequency, waveguide, Mfcn)
    LMP.solvemodalequation(ea, modeequation)
end


@testset "modefinderX.jl" begin
    @info "Testing modefinderX"

    @test test_dX()
    @test test_R2XX2R()
    @test test_modalequationX()
end
