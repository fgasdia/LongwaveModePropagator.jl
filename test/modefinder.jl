function test_physicalmodeequation(scenario)
    @unpack ea, tx, bfield, species, ground = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)

    me = LMP.PhysicalModeEquation(ea, tx.frequency, waveguide)
    @test me isa LMP.ModeEquation
    @test me isa LMP.PhysicalModeEquation

    me2 = LMP.PhysicalModeEquation(tx.frequency, waveguide)
    @test me2 isa LMP.PhysicalModeEquation
    @test iszero(me2.ea.θ)

    @test LMP.setea(ea, me2) == me
    @test LMP.setea(ea.θ, me2) == me
end

function test_roots(scenario)
    x = 0.005
    z = complex(x, -x)

    @test LMP.isroot(x)
    @test LMP.isroot(x, atol=1e-6) == false
    @test LMP.isroot(z)
    @test LMP.isroot(z, atol=1e-6) == false

    z2 = complex(1, 0)
    z3 = complex(0, 1)
    @test LMP.isroot(z2) == false
    @test LMP.isroot(z3) == false

    @unpack ea, tx, bfield, species, ground = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    me = LMP.PhysicalModeEquation(tx.frequency, waveguide)

    roots = copy(TEST_ROOTS)
    @test LMP.filterroots!(roots, me) == roots
    @test LMP.filterroots!(roots, tx.frequency, waveguide) == roots

    push!(roots, complex(1.5, -0.5))
    @test LMP.filterroots!(roots, me) == TEST_ROOTS
    @test LMP.filterroots!(roots, tx.frequency, waveguide) == TEST_ROOTS

    roots2 = copy(roots)
    roots2[1] += 0.001
    f = LMP.solvemodalequation(LMP.setea(roots2[1], me))
    @test LMP.isroot(f, atol=0.5)  # NOTE: atol is for value of modal equation, not θ
    @test LMP.filterroots!(roots2, me, atol=0.5) == roots2
    @test LMP.filterroots!(roots2, me) == roots[2:end]
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
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    e = LMP.bookerwavefields(ea, M)
    @test LMP.bookerreflection(ea, e) == LMP.bookerreflection(ea, M)

    #==
    Iterative solution for sharp boundary reflection matrix
    ==#
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
        @test err_func(getindex.(initRs, n), getindex.(Rs, n)) < 1e-8
    end
end

function test_bookerreflection_vertical(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)
    R = LMP.bookerreflection(ea, M)

    @test R[1,2] ≈ R[2,1]
end

function test_dbookerreflection(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    # Iterative solution
    initdRs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(θs))
    dRs = similar(initdRs)
    for i in eachindex(θs)
        ea = EigenAngle(θs[i])
        T = LMP.tmatrix(ea, M)
        W = LMP.wmatrix(ea, T)
        dW = LMP.wmatrix(ea, M, T, LMP.Dθ())

        R = LMP.bookerreflection(ea, M)

        initdR = LMP.bookerreflection(ea, M, LMP.Dθ())[2]
        resdR = nlsolve((f,dR)->sharpR!(f,dR,R,dW,W), Array(initdR))
        dR = SMatrix{2,2}(resdR.zero)

        initdRs[i] = initdR
        dRs[i] = dR
    end

    for n in 1:4
        @test err_func(getindex.(initdRs, n), getindex.(dRs, n)) < 1e-6  # TODO BUG I suspect a bug because iterative matches perfectly
    end

    #==
    Finit difference derivative
    ==#
    for i = 1:4
        Rref(θ) = (ea = EigenAngle(θ); LMP.bookerreflection(ea, M)[i])
        dRref = FiniteDiff.finite_difference_derivative(Rref, θs, Val{:central})
        R(θ) = (ea = EigenAngle(θ); LMP.bookerreflection(ea, M, LMP.Dθ())[1][i])
        dR(θ) = (ea = EigenAngle(θ); LMP.bookerreflection(ea, M, LMP.Dθ())[2][i])

        @test Rref.(θs) ≈ R.(θs)
        @test err_func(dRref, dR.(θs)) < 1e-6
    end
end

function wmatrix_deriv(scenario)
    @unpack tx, bfield, species = scenario

    M = LMP.susceptibility(70e3, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Wfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M); LMP.wmatrix(ea, T)[i][j])
            dWref = FiniteDiff.finite_difference_derivative(Wfcn, θs, Val{:central})
            dW(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M); LMP.wmatrix(ea, M, T, LMP.Dθ())[i][j])

            err_func(dW.(θs), dWref) < 1e-3 || return false
        end
    end

    return true
end

function integratedreflection_deriv(scenario)
    @unpack tx, ground, bfield, species = scenario
    freq = tx.frequency

    # params = LMPParams(integrationparams=IntegrationParams(1e-7, LMP.RK4(), false))
    params = LMPParams()
    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(freq, waveguide)

    Rref(θ,m) = (m = LMP.setea(EigenAngle(θ), m); LMP.integratedreflection(m, params=params))
    RdR(θ,m) = (m = LMP.setea(EigenAngle(θ), m); LMP.integratedreflection(m, LMP.Dθ(), params=params))

    Rs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(θs))
    dRs = similar(Rs)
    Rrefs = similar(Rs)
    dRrefs = similar(Rs)
    Threads.@threads for i in 1:length(θs)  # NOTE: this also effectively checks for thread safety
        v = RdR(θs[i], modeequation)
        Rs[i] = v[SVector(1,2),:]
        dRs[i] = v[SVector(3,4),:]
        Rrefs[i] = Rref(θs[i], modeequation)
        dRrefs[i] = FiniteDiff.finite_difference_derivative(z->Rref(z, modeequation), θs[i],
                                                            Val{:central})
    end

    for i = 1:4
        R = [v[i] for v in Rs[1:end-1]]
        dR = [v[i] for v in dRs[1:end-1]]
        Rr = [v[i] for v in Rrefs[1:end-1]]
        dRr = [v[i] for v in dRrefs[1:end-1]]

        isapprox(R, Rr, rtol=1e-5) || return false
        isapprox(dR, dRr, rtol=1e-3) || return false
    end

    return true
end

function fresnelreflection_deriv(scenario)
    @unpack tx, ground = scenario
    freq = tx.frequency

    for i = 1:4
        Rgref(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq)[i])
        dRgref = FiniteDiff.finite_difference_derivative(Rgref, θs, Val{:central})
        Rg(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq, LMP.Dθ())[1][i])
        dRg(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq, LMP.Dθ())[2][i])

        Rgref.(θs) ≈ Rg.(θs) || return false
        err_func(dRg.(θs), dRgref) < 1e-6 || return false
    end

    return true
end

function modalequation_deriv(scenario)
    @unpack ea, tx, ground, bfield, species = scenario
    freq = tx.frequency

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(ea, freq, waveguide)

    dFref = FiniteDiff.finite_difference_derivative(z->LMP.solvemodalequation(z, modeequation),
        θs, Val{:central})

    dFs = Vector{ComplexF64}(undef, length(θs))
    Threads.@threads for i in eachindex(θs)
        dFdθ, R, Rg = LMP.solvemodalequation(θs[i], modeequation, LMP.Dθ())
        dFs[i] = dFdθ
    end

    return isapprox(dFs, dFref, rtol=1e-4)
end

########
# Non-derivative
########

function test_wmatrix(scenario)
    # Check that "analytical" solution for W matches numerical
    @unpack ea, tx, bfield, species = scenario

    C = ea.cosθ
    L = [C 0 -C 0;
         0 -1 0 -1;
         0 -C 0 C;
         1 0 1 0]

    M = LMP.susceptibility(80e3, tx.frequency, bfield, species)
    T = LMP.tmatrix(ea, M)
    W = LMP.wmatrix(ea, T)

    return [W[1] W[3]; W[2] W[4]] ≈ 2*(L\T)*L
end



function verticalreflection(scenario)
    @unpack ea, tx, ground, bfield, species = scenario

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(tx.frequency, waveguide)

    R = LMP.integratedreflection(modeequation)

    return R[1,2] ≈ R[2,1]
end

function pecground()
    pec_ground = LMP.Ground(1, 1e12)
    vertical_ea = LMP.EigenAngle(π/2)

    Rg = LMP.fresnelreflection(vertical_ea, pec_ground, Frequency(24e3))

    return isapprox(abs.(Rg), I, atol=1e-7)
end

function modalequation(scenario)
    @unpack ea, tx, ground, bfield, species = scenario

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(ea, tx.frequency, waveguide)

    f = LMP.solvemodalequation(modeequation)

    return LMP.isroot(f)
end

function modefinder(scenario)
    @unpack tx, bfield, species, ground = scenario
    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(tx.frequency, waveguide)

    origcoords = LMP.defaultcoordinates(tx.frequency)

    params = LMPParams(grpfparams=LMP.GRPFParams(100000, 1e-6, true))
    modes = LMP.findmodes(modeequation, origcoords, params=params)

    for m in modes
        f = LMP.solvemodalequation(m, modeequation, params=params)
        LMP.isroot(f) || return false
    end
    # return modes
    return true
end

function evalroot(root, scenario)
    @unpack tx, bfield, species, ground = scenario
    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)

    modeequation = LMP.PhysicalModeEquation(tx.frequency, waveguide)
    LMP.solvemodalequation(EigenAngle(root), modeequation)
end


@testset "modefinder.jl" begin
    @info "Testing modefinder"

    @test pecground()

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        test_physicalmodeequation(scn)  # just test with one scenario?

        test_bookerreflection(scn)
        test_dbookerreflection(scn)

        @test test_wmatrix(scn)
        @test modefinder(scn)
        @test iterativesharpboundary(scn)
    end
    for scn in (resonant_scenario, )
        test_roots(scn)

        @test modalequation(scn)
    end
    for scn in (verticalB_scenario, )
        test_bookerreflection_vertical(scn)

        @test verticalreflection(scn)
    end

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        @test wmatrix_deriv(scn)
        @test fresnelreflection_deriv(scn)
        @test integratedreflection_deriv(scn) # again, a problem with off-diag
        @test modalequation_deriv(scn)
    end
end
