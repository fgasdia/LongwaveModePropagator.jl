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
    susceptibilityfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    Mtop = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)
    Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
    Xtop = LWMS.R2X(ea, Rtop)

    dzparams = LWMS.DZParams(ea, tx.frequency, susceptibilityfcn)
    X = LWMS.dXdz(Xtop, dzparams, 80e3)
    Xref = dXdzW(Xtop, dzparams, 80e3)

    return X ≈ Xref
end

function test_R2XX2R()
    ea = EigenAngle(deg2rad(85.0-1.0im))
    R = rand(SMatrix{2,2,ComplexF64})
    X = LWMS.R2X(ea, R)
    Rtest = LWMS.X2R(ea, X)

    return Rtest ≈ X
end

function test_modefinderX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    susceptibilityfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = GRPFParams(est_num_nodes, 1e-8, true)

    @unpack bfield, species = waveguide
    susceptibilityfcn = LWMS.susceptibilityinterpolator(tx.frequency, bfield, electrons)

    roots, poles = grpf(z->solvemodalequationX(EigenAngle(z), tx.frequency, waveguide, susceptibilityfcn),
                        origcoords, grpfparams)

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
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    susceptibilityfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = GRPFParams(est_num_nodes, 1e-8, true)

    @unpack bfield, species = waveguide
    # susceptibilityfcn = LWMS.susceptibilityinterpolator(tx.frequency, bfield, species)
    susceptibilitymfcn = alt->LWMS.susceptibility(alt, tx.frequency, bfield, species)

    LWMS.integratedreflectionX(ea, tx.frequency, waveguide, susceptibilitymfcn)
end

function test_modalequationX()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    susceptibilityfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    LWMS.solvemodalequationX(ea, tx.frequency, waveguide, susceptibilityfcn)

    # eas = EigenAngle.(rand(ComplexF64, 1000))
    # for i in eachindex(eas)
    #     LWMS.solvemodalequationX(eas[i], tx.frequency, waveguide, susceptibilityfcn)
    # end
end

function modefinderX()
    ea, tx, ground, bfield, electrons = scenario()

    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    # origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    # est_num_nodes = ceil(Int, length(origcoords)*1.5)
    # grpfparams = GRPFParams(est_num_nodes, 1e-8, true)

    Zb = deg2rad(complex(40.0, -10.0))
    Ze = deg2rad(complex(89.9, 0.0))
    d = deg2rad(60)
    Δr = deg2rad(0.6)
    origcoords = LWMS.uppertriangularrectdomain(Zb, Ze, Δr)

    # Δr from 0.5->0.25 => time from 3.8->5.3 sec
    # tolerance from 1e-8->1e-7 => time from 5.3->4.6 sec
    # origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    # est_num_nodes = ceil(Int, length(origcoords)*1.5)
    # grpfparams = GRPFParams(est_num_nodes, 1e-7, false)

    # if approximate
    #     Mfcn = susceptibilityinterpolator(frequency, bfield, species)
    # else
        # Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, electrons)
    # end

    # for i in eachindex(origcoords)
    #     LWMS.solvemodalequation(EigenAngle(origcoords[i]), tx.frequency, waveguide, Mfcn)
    # end

    # f(z) = LWMS.solvemodalequationX(EigenAngle(z), tx.frequency, waveguide, Mfcn)
    # f(80e3)
    # roots, poles = grpf(f, origcoords, grpfparams)

    # return roots
    modes = LWMS.findmodes(origcoords, tx.frequency, waveguide, grpfparams, true, true)

    M(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, electrons)
    for m in modes
        fval = LWMS.solvemodalequationX(m, tx.frequency, waveguide, M)
        # println(f)
        LWMS.isroot(fval) || return false
    end
    return true
end


function tdx()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    susceptibilityfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    dzparams = LWMS.DZParams(ea, tx.frequency, susceptibilityfcn)

    Xtop = rand(SMatrix{2,2,ComplexF64})
    # Xtop = @SMatrix [2.0+1.0im 2.1+1.1im
                     # 19+0.2im 13.0+0.5im]
    X = LWMS.dXdz(Xtop, dzparams, 80e3)
    # R = LWMS.dRdz(Xtop, dzparams, 80e3)
end

function tsx()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    susceptibilityfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    # dzparams = LWMS.DZParams(ea, tx.frequency, susceptibilityfcn)

    LWMS.solvemodalequationX(ea, tx.frequency, waveguide, susceptibilityfcn)
    # LWMS.solvemodalequation(ea, tx.frequency, waveguide, susceptibilityfcn)
end


@testset "modefinderX.jl" begin
    @info "Testing modefinderX"

    @test test_dX()
    @test test_R2XX2R()
end
