function test_waitsusceptibilityinterp(scenario)
    @unpack tx, bfield, species = scenario

    zs = rand(100000).*(LWMS.TOPHEIGHT-LWMS.BOTTOMHEIGHT) .+ LWMS.BOTTOMHEIGHT
    trueM = LWMS.susceptibility.(zs, (tx.frequency,), (bfield,), (species,))

    itp = LWMS.susceptibilityinterpolator(tx.frequency, bfield, species)
    interpM = itp(zs)

    # In use, recommend
    # import FunctionWrappers: FunctionWrapper
    # itp = LWMS.susceptibilityinterpolator(frequency, bfield, species)
    # Mfcn = FunctionWrapper{SMatrix{3,3,ComplexF64,9}, Tuple{Float64}}(itp)

    return trueM ≈ interpM
end

function evalMfcn(Mfcn, zs)
    for i in eachindex(zs)
        Mfcn(zs[i])
    end
end

function test_bookerquarticM(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, species)
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

function test_bookerquarticT(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, species)
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

function test_bookerquartics(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    qM, BM = LWMS.bookerquartic(ea, M)
    qT, BT = LWMS.bookerquartic(T)

    return qM ≈ qT
end

@testset "magnetoionic.jl" begin
    @info "Testing magnetoionic"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        @test test_waitsusceptibilityinterp(scn)
        @test test_waitsusceptibilityinterp(scn)
        @test_skip test_nonwaitsusceptibilityinterp(scn)

        # TODO: Off-diagonal terms should be 0 with no B field
        @test_skip +(M[1,2], M[1,3], M[2,1], M[2,3], M[3,1], M[3,2]) == 0

        @test test_bookerquarticM(scn)
        @test test_bookerquarticT(scn)
        @test test_bookerquartics(scn)
    end
end
