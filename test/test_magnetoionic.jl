function scenario()
    ea = EigenAngle(deg2rad(complex(85.0, -1.0)))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    # `tx` isn't just bits
    ground = Ground(15, 0.001)
    bfield = BField(50e-6, deg2rad(90), 0)

    electrons = Species(QE, ME,
                        z -> waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
                        z -> electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT))

    return ea, tx, ground, bfield, electrons
end

function test_waitsusceptibilityinterp()
    ea, tx, ground, bfield, electrons = scenario()

    zs = rand(100000).*(LWMS.TOPHEIGHT-LWMS.BOTTOMHEIGHT) .+ LWMS.BOTTOMHEIGHT
    trueM = LWMS.susceptibility.(zs, (tx.frequency,), (bfield,), (electrons,))

    itp = LWMS.susceptibilityinterpolator(tx.frequency, bfield, electrons)
    interpM = itp(zs)

    return trueM ≈ interpM
end

function evalMfcn(Mfcn, zs)
    for i in eachindex(zs)
        Mfcn(zs[i])
    end
end

function test_bookerquarticM()
    ea, tx, ground, bfield, electrons = scenario()

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)
    q, B = LWMS.bookerquartic(ea, M)

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1] M[1,2] S*q[i]+M[1,3];
             M[2,1] 1-q[i]^2-S²+M[2,2] M[2,3];
             S*q[i]+M[3,1] M[3,2] C²+M[3,3]]
        LWMS.isroot(det(G), atol=sqrt(eps())) || return false
    end

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        LWMS.isroot(booker, atol=sqrt(eps())) || return false
    end

    return true
end

function test_bookerquarticT()
    ea, tx, ground, bfield, electrons = scenario()

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)
    q, B = LWMS.bookerquartic(T)

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1] M[1,2] S*q[i]+M[1,3];
             M[2,1] 1-q[i]^2-S²+M[2,2] M[2,3];
             S*q[i]+M[3,1] M[3,2] C²+M[3,3]]
        LWMS.isroot(det(G), atol=sqrt(eps())) || return false
    end

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        LWMS.isroot(booker, atol=sqrt(eps())) || return false
    end

    return true
end

function test_bookerquartics()
    ea, tx, ground, bfield, electrons = scenario()

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)

    qM, BM = LWMS.bookerquartic(ea, M)
    qT, BT = LWMS.bookerquartic(T)

    return qM ≈ qT
end

@testset "magnetoionic.jl" begin
    @info "Testing magnetoionic"

    @test test_waitsusceptibilityinterp()
    @test_skip test_nonwaitsusceptibilityinterp()

    # TODO: Off-diagonal terms should be 0 with no B field
    @test_skip +(M[1,2], M[1,3], M[2,1], M[2,3], M[3,1], M[3,2]) == 0

    @test test_bookerquarticM()
    @test test_bookerquarticT()
    @test test_bookerquartics()
end
