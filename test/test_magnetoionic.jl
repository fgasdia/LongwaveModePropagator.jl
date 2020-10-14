function test_bookerquarticM(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    q, B = LWMS.bookerquartic(ea, M)

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1]    M[1,2]              S*q[i]+M[1,3];
             M[2,1]             1-q[i]^2-S²+M[2,2]  M[2,3];
             S*q[i]+M[3,1]      M[3,2]              C²+M[3,3]]

        LWMS.isroot(det(G), atol=1e-5) || return false
    end

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        LWMS.isroot(booker, atol=1e-5) || return false
    end

    return true
end

function test_bookerquarticT(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)
    q, B = LWMS.bookerquartic(T)

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1]    M[1,2]              S*q[i]+M[1,3];
             M[2,1]             1-q[i]^2-S²+M[2,2]  M[2,3];
             S*q[i]+M[3,1]      M[3,2]              C²+M[3,3]]
        LWMS.isroot(det(G), atol=1e-5) || return false
    end

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        LWMS.isroot(booker, atol=1e-5) || return false
    end

    return true
end

function test_bookerquartics(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    qM, BM = LWMS.bookerquartic(ea, M)
    qT, BT = LWMS.bookerquartic(T)

    return qM ≈ qT
end

function test_bookerquarticM_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)

    for i = 1:4
        qfcn(θ) = (ea = EigenAngle(θ); (q, B) = LWMS.bookerquartic(ea, M);
            sort!(q, by=LWMS.upgoing); q[i])
        dqref = FiniteDiff.finite_difference_derivative(qfcn, θs, Val{:central})
        dq(θ) = (ea = EigenAngle(θ); (q, B) = LWMS.bookerquartic(ea, M);
            sort!(q, by=LWMS.upgoing); LWMS.bookerquartic(ea, M, q, B, LWMS.Dθ())[i])

        err_func(dq.(θs), dqref) < 1e-6 || return false
    end

    return true
end

function tmatrix_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Tfcn(θ) = (ea = EigenAngle(θ); T = LWMS.tmatrix(ea, M)[i,j])
            dTref = FiniteDiff.finite_difference_derivative(Tfcn, θs, Val{:central})
            dT(θ) = (ea = EigenAngle(θ); T = LWMS.tmatrix(ea, M, LWMS.Dθ())[i,j])

            err_func(dT.(θs), dTref) < 1e-6 || return false
        end
    end

    return true
end

@testset "magnetoionic.jl" begin
    @info "Testing magnetoionic"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        @test tmatrix_deriv(scn)

        # TODO: Off-diagonal terms should be 0 with no B field
        @test_skip +(M[1,2], M[1,3], M[2,1], M[2,3], M[3,1], M[3,2]) == 0

        @test test_bookerquarticM(scn)
        @test test_bookerquarticT(scn)
        @test test_bookerquartics(scn)
        @test test_bookerquarticM_deriv(scn)
    end
end
