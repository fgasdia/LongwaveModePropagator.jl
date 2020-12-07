function test_susceptibility(scenario)
    @unpack ea, tx, bfield, species, ground = scenario()

    M1 = LMP.susceptibility(70e3, tx.frequency, bfield, species)
    M2 = LMP.susceptibility(70e3, tx.frequency, bfield, species, params=LMPParams())
    M3 = LMP.susceptibility(70e3, tx.frequency, bfield, species,
                             params=LMPParams(earthradius=6350e3))
    @test M1 == M2
    @test !(M2 ≈ M3)

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(tx.frequency, waveguide)

    M4 = LMP.susceptibility(70e3, tx.frequency, waveguide)
    M5 = LMP.susceptibility(70e3, tx.frequency, waveguide, params=LMPParams())
    M6 = LMP.susceptibility(70e3, tx.frequency, waveguide,
                             params=LMPParams(earthradius=6350e3))

    M7 = LMP.susceptibility(70e3, modeequation)
    M8 = LMP.susceptibility(70e3, modeequation, params=LMPParams())
    M9 = LMP.susceptibility(70e3, modeequation,
                             params=LMPParams(earthradius=6350e3))

    @test M4 == M5 == M1
    @test M6 == M3

    @test M7 == M8 == M1
    @test M9 == M3
end

function tmatrix_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Tfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M)[i,j])
            dTref = FiniteDiff.finite_difference_derivative(Tfcn, θs, Val{:central})
            dT(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M, LMP.Dθ())[i,j])

            @test err_func(dT.(θs), dTref) < 1e-6
        end
    end
end

function test_bookerquarticM(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)
    q, B = LMP.bookerquartic(ea, M)

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1]    M[1,2]              S*q[i]+M[1,3];
             M[2,1]             1-q[i]^2-S²+M[2,2]  M[2,3];
             S*q[i]+M[3,1]      M[3,2]              C²+M[3,3]]

        @test LMP.isroot(det(G), atol=1e-5)  # real barely fails 1e-6 for resonant_scenario
    end

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        @test LMP.isroot(booker, atol=1e-6)
    end
end

function test_bookerquarticT(scenario)
    @unpack ea, tx, bfield, species = scenario()

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
        @test LMP.isroot(det(G), atol=1e-6)
    end

    # eigvals is >20 times slower than bookerquartic
    @test sort(eigvals(Array(T)), by=LMP.upgoing) ≈ sort(q, by=LMP.upgoing)

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        @test LMP.isroot(booker, atol=1e-6)
    end
end

function test_bookerquartics(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)
    T = LMP.tmatrix(ea, M)

    qM, BM = LMP.bookerquartic(ea, M)
    qT, BT = LMP.bookerquartic(T)

    sort!(qM, by=LMP.upgoing)
    sort!(qT, by=LMP.upgoing)

    @test qM ≈ qT
end

function test_bookerquarticM_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    for i = 1:4
        qfcn(θ) = (ea = EigenAngle(θ); (q, B) = LMP.bookerquartic(ea, M);
            sort!(q, by=LMP.upgoing); q[i])
        dqref = FiniteDiff.finite_difference_derivative(qfcn, θs, Val{:central})
        dq(θ) = (ea = EigenAngle(θ); (q, B) = LMP.bookerquartic(ea, M);
            sort!(q, by=LMP.upgoing); LMP.dbookerquartic(ea, M, q, B)[i])

        @test err_func(dq.(θs), dqref) < 1e-6
    end
end

function test_bookerquarticT_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    for i = 1:4
        qfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M);
                   (q, B) = LMP.bookerquartic(T); LMP.sortquarticroots!(q); q[i])
        dqref = FiniteDiff.finite_difference_derivative(qfcn, θs, Val{:central})
        dq(θ) = (ea = EigenAngle(θ);
                 T = LMP.tmatrix(ea, M); dT = LMP.tmatrix(ea, M, LMP.Dθ());
                 (q, B) = LMP.bookerquartic(T);
                 LMP.sortquarticroots!(q); LMP.dbookerquartic(T, dT, q, B)[i])

        @test err_func(dq.(θs), dqref) < 1e-6
    end
end

@testset "magnetoionic.jl" begin
    @info "Testing magnetoionic"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        test_susceptibility(scn)

        tmatrix_deriv(scn)

        test_bookerquarticM(scn)
        test_bookerquarticT(scn)
        test_bookerquartics(scn)
        test_bookerquarticM_deriv(scn)
        test_bookerquarticT_deriv(scn)
    end
end
