function test_bookerquarticM(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)
    q, B = LMP.bookerquartic(ea, M)

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1]    M[1,2]              S*q[i]+M[1,3];
             M[2,1]             1-q[i]^2-S²+M[2,2]  M[2,3];
             S*q[i]+M[3,1]      M[3,2]              C²+M[3,3]]

        @test isroot(det(G); atol=1e-5)  # real barely fails 1e-6 for resonant_scenario
    end

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        @test isroot(booker; atol=1e-6)
    end
end

function test_bookerquarticT(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)
    T = LMP.tmatrix(ea, M)
    q, B = LMP.bookerquartic(T)

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1]    M[1,2]              S*q[i]+M[1,3];
             M[2,1]             1-q[i]^2-S²+M[2,2]  M[2,3];
             S*q[i]+M[3,1]      M[3,2]              C²+M[3,3]]
        @test isroot(det(G); atol=1e-6)
    end

    # eigvals is >20 times slower than bookerquartic
    @test sort(eigvals(Array(T)); by=LMP.upgoing) ≈ sort(q; by=LMP.upgoing)

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        @test isroot(booker; atol=1e-6)
    end
end

function test_bookerquartics(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)
    T = LMP.tmatrix(ea, M)

    qM, BM = LMP.bookerquartic(ea, M)
    qT, BT = LMP.bookerquartic(T)

    sort!(qM; by=LMP.upgoing)
    sort!(qT; by=LMP.upgoing)

    @test qM ≈ qT
end

function test_bookerquarticM_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    for i = 1:4
        qfcn(θ) = (ea = EigenAngle(θ); (q, B) = LMP.bookerquartic(ea, M);
            sort!(q; by=LMP.upgoing); q[i])
        dqref = FiniteDiff.finite_difference_derivative(qfcn, θs, Val{:central})
        dq(θ) = (ea = EigenAngle(θ); (q, B) = LMP.bookerquartic(ea, M);
            sort!(q; by=LMP.upgoing); LMP.dbookerquartic(ea, M, q, B)[i])

        @test maxabsdiff(dq.(θs), dqref) < 1e-6
    end
end

function test_bookerquarticT_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    for i = 1:4
        qfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M);
                   (q, B) = LMP.bookerquartic(T); LMP.sortquarticroots!(q); q[i])
        dqref = FiniteDiff.finite_difference_derivative(qfcn, θs, Val{:central})
        dq(θ) = (ea = EigenAngle(θ);
                 T = LMP.tmatrix(ea, M); dT = LMP.dtmatrix(ea, M);
                 (q, B) = LMP.bookerquartic(T);
                 LMP.sortquarticroots!(q); LMP.dbookerquartic(T, dT, q, B)[i])

        @test maxabsdiff(dq.(θs), dqref) < 1e-6
    end
end

function test_bookerwavefields(scenario)
    @unpack ea, bfield, tx, species = scenario

    topheight = first(LMPParams().wavefieldheights)
    M = LMP.susceptibility(topheight, tx.frequency, bfield, species)
    T = LMP.tmatrix(ea, M)

    # Verify bookerwavefields produces a valid solution
    e = LMP.bookerwavefields(T)
    q, B = LMP.bookerquartic(T)
    LMP.sortquarticroots!(q)

    for i = 1:2
        @test T*e[:,i] ≈ q[i]*e[:,i]
    end

    e2 = LMP.bookerwavefields(ea, M)
    @test e ≈ e2
end

function test_bookerwavefields_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    # Compare analytical solution with finite diff
    for i = 1:4
        eref(θ) = (ea = EigenAngle(θ); LMP.bookerwavefields(ea, M)[i])
        deref = FiniteDiff.finite_difference_derivative(eref, θs, Val{:central})
        e(θ) = (ea = EigenAngle(θ); LMP.bookerwavefields(ea, M, LMP.Dθ())[1][i])
        de(θ) = (ea = EigenAngle(θ); LMP.bookerwavefields(ea, M, LMP.Dθ())[2][i])

        @test eref.(θs) ≈ e.(θs)
        @test maxabsdiff(deref, de.(θs)) < 1e-3
    end

    # Compare functions using `M` and `T`
    de1 = Vector{SMatrix{4,2,ComplexF64,8}}(undef, length(θs))
    de2 = similar(de1)
    for i in eachindex(θs)
        ea = EigenAngle(θs[i])
        T = LMP.tmatrix(ea, M)
        dT = LMP.dtmatrix(ea, M)
        de1[i] = LMP.bookerwavefields(ea, M, LMP.Dθ())[2]
        de2[i] = LMP.bookerwavefields(T, dT, LMP.Dθ())[2]
    end
    @test de1 ≈ de2
end

function test_bookerreflection_vertical(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)
    R = LMP.bookerreflection(ea, M)

    @test R[1,2] ≈ R[2,1]
end

function sharpR!(f, R, W)
    #==
    The right side of ``2i R′ = w₂₁ + w₂₂R - Rw₁₁ - Rw₁₂R``, which can be used to solve for
    reflection from the sharply bounded ionosphere where ``dR/dz = 0``.
    ==#
    # See [^Budden1988] sec. 18.10.
    f .= W[2] + W[4]*R - R*W[1] - R*W[3]*R
end

function sharpR!(f, dR, R, dW, W)
    # ``dR/dz/dθ`` which can be used to solve for ``dR/dθ`` from a sharply bounded ionosphere.
    f .= dW[2] + W[4]*dR + dW[4]*R - R*dW[1] - dR*W[1] - R*W[3]*dR - R*dW[3]*R - dR*W[3]*R
end

function test_bookerreflection(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    e = LMP.bookerwavefields(ea, M)
    @test LMP.bookerreflection(ea, e) == LMP.bookerreflection(ea, M)

    # Iterative solution for sharp boundary reflection matrix
    initRs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(θs))
    Rs = similar(initRs)
    for i in eachindex(θs)
        ea = EigenAngle(θs[i])
        T = LMP.tmatrix(ea, M)
        W = LMP.wmatrix(ea, T)

        initR = LMP.bookerreflection(ea, M)
        res = nlsolve((f,R)->sharpR!(f,R,W), Array(initR))
        R = SMatrix{2,2}(res.zero)

        initRs[i] = initR
        Rs[i] = R
    end

    for n in 1:4
        @test maxabsdiff(getindex.(initRs, n), getindex.(Rs, n)) < 1e-8
    end

    # Matrix solution
    e = LMP.bookerwavefields(ea, M)
    wavefieldR = LMP.bookerreflection(ea, e)

    Cinv = ea.secθ
    Sv_inv = SMatrix{4,4}(Cinv, 0, -Cinv, 0,
                          0, -1, 0, -1,
                          0, -Cinv, 0, Cinv,
                          1, 0, 1, 0)

    f1 = Sv_inv*e[:,1]
    f2 = Sv_inv*e[:,2]

    U = SMatrix{2,2}(f1[1], f1[2], f2[1], f2[2])
    D = SMatrix{2,2}(f1[3], f1[4], f2[3], f2[4])

    @test D/U ≈ wavefieldR
end

function test_dbookerreflection(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    # Iterative solution
    initdRs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(θs))
    dRs = similar(initdRs)
    for i in eachindex(θs)
        ea = EigenAngle(θs[i])
        T = LMP.tmatrix(ea, M)
        dT = LMP.dtmatrix(ea, M)
        W = LMP.wmatrix(ea, T)
        dW = LMP.dwmatrix(ea, T, dT)

        R = LMP.bookerreflection(ea, M)

        initdR = LMP.bookerreflection(ea, M, LMP.Dθ())[2]
        resdR = nlsolve((f,dR)->sharpR!(f,dR,R,dW,W), Array(initdR))
        dR = SMatrix{2,2}(resdR.zero)

        initdRs[i] = initdR
        dRs[i] = dR
    end

    for n in 1:4
        @test maxabsdiff(getindex.(initdRs, n), getindex.(dRs, n)) < 1e-4
    end

    # Finite difference derivative
    for i = 1:4
        Rref(θ) = (ea = EigenAngle(θ); LMP.bookerreflection(ea, M)[i])
        dRref = FiniteDiff.finite_difference_derivative(Rref, θs, Val{:central})
        R(θ) = (ea = EigenAngle(θ); LMP.bookerreflection(ea, M, LMP.Dθ())[1][i])
        dR(θ) = (ea = EigenAngle(θ); LMP.bookerreflection(ea, M, LMP.Dθ())[2][i])

        @test Rref.(θs) ≈ R.(θs)
        @test maxabsdiff(dRref, dR.(θs)) < 1e-4
    end
end

@testset "bookerquartic.jl" begin
    @info "Testing bookerquartic"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        multiplespecies_scenario)
        test_bookerquarticM(scn)
        test_bookerquarticT(scn)
        test_bookerquartics(scn)
        test_bookerquarticM_deriv(scn)
        test_bookerquarticT_deriv(scn)

        test_bookerwavefields(scn)
        test_bookerwavefields_deriv(scn)

        test_bookerreflection(scn)
        test_dbookerreflection(scn)
    end

    for scn in (verticalB_scenario, )
        test_bookerreflection_vertical(scn)
    end
end
