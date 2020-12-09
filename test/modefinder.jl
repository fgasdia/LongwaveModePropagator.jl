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

    # isroot Real and Complex
    @test LMP.isroot(x)
    @test LMP.isroot(x, atol=1e-6) == false
    @test LMP.isroot(z)
    @test LMP.isroot(z, atol=1e-6) == false

    z2 = complex(1, 0)
    z3 = complex(0, 1)
    @test LMP.isroot(z2) == false
    @test LMP.isroot(z3) == false

    # filterroots! setup
    @unpack ea, tx, bfield, species, ground = scenario()
    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    me = LMP.PhysicalModeEquation(tx.frequency, waveguide)

    # nothing filtered
    roots = copy(TEST_ROOTS)
    @test LMP.filterroots!(roots, me) == roots
    @test LMP.filterroots!(roots, tx.frequency, waveguide) == roots

    # filter bad root
    push!(roots, complex(1.5, -0.5))
    @test LMP.filterroots!(roots, me) == TEST_ROOTS
    @test LMP.filterroots!(roots, tx.frequency, waveguide) == TEST_ROOTS

    # filter (or not) with different tolerance
    roots2 = copy(roots)
    roots2[1] += 0.001
    f = LMP.solvemodalequation(LMP.setea(roots2[1], me))
    @test LMP.isroot(f, atol=0.5)  # NOTE: atol is for value of modal equation, not θ
    @test LMP.filterroots!(roots2, me, atol=0.5) == roots2
    @test LMP.filterroots!(roots2, me) == roots[2:end]
end

function test_wmatrix(scenario)
    @unpack ea, tx, bfield, species = scenario()

    C = ea.cosθ
    L = @SMatrix [C 0 -C 0;
                  0 -1 0 -1;
                  0 -C 0 C;
                  1 0 1 0]

    M = LMP.susceptibility(80e3, tx.frequency, bfield, species)
    T = LMP.tmatrix(ea, M)
    W = LMP.wmatrix(ea, T)

    @test [W[1] W[3]; W[2] W[4]] ≈ 2*(L\T)*L
end

function test_wmatrix_deriv(scenario)
    @unpack tx, bfield, species = scenario()

    M = LMP.susceptibility(80e3, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Wfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M); LMP.wmatrix(ea, T)[i][j])
            dWref = FiniteDiff.finite_difference_derivative(Wfcn, θs, Val{:central})
            dW(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M);
                     dT = LMP.dtmatrix(ea, M); LMP.dwmatrix(ea, T, dT)[i][j])

            @test err_func(dW.(θs), dWref) < 1e-3
        end
    end
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

function test_fresnelreflection(scenario)
    @unpack ea, tx, bfield, species, ground = scenario()

    # PEC ground
    pec_ground = LMP.Ground(1, 1e12)
    vertical_ea = LMP.EigenAngle(π/2)
    Rg = LMP.fresnelreflection(vertical_ea, pec_ground, Frequency(24e3))
    @test isapprox(abs.(Rg), I, atol=1e-7)

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    me = LMP.PhysicalModeEquation(ea, tx.frequency, waveguide)
    @test LMP.fresnelreflection(ea, ground, tx.frequency) == LMP.fresnelreflection(me)
end

function test_fresnelreflection_deriv(scenario)
    @unpack ea, tx, bfield, species, ground = scenario()
    freq = tx.frequency

    for i = 1:4
        Rgref(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq)[i])
        dRgref = FiniteDiff.finite_difference_derivative(Rgref, θs, Val{:central})
        Rg(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq, LMP.Dθ())[1][i])
        dRg(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq, LMP.Dθ())[2][i])

        @test Rgref.(θs) ≈ Rg.(θs)
        @test err_func(dRg.(θs), dRgref) < 1e-6
    end

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    me = LMP.PhysicalModeEquation(ea, tx.frequency, waveguide)
    Rg1 = LMP.fresnelreflection(ea, ground, freq, LMP.Dθ())
    Rg2 = LMP.fresnelreflection(me, LMP.Dθ())
    @test Rg1 == Rg2
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

function verticalreflection(scenario)
    @unpack ea, tx, ground, bfield, species = scenario

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(tx.frequency, waveguide)

    R = LMP.integratedreflection(modeequation)

    return R[1,2] ≈ R[2,1]
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

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        test_physicalmodeequation(scn)  # just test with one scenario?

        test_wmatrix(scn)
        test_wmatrix_deriv(scn)

        test_fresnelreflection(scn)
        test_fresnelreflection_deriv(scn)

        @test modefinder(scn)
        @test iterativesharpboundary(scn)
        @test integratedreflection_deriv(scn) # again, a problem with off-diag
        @test modalequation_deriv(scn)
    end
    for scn in (resonant_scenario, )
        test_roots(scn)

        @test modalequation(scn)
    end
    for scn in (verticalB_scenario, )
        @test verticalreflection(scn)
    end
end
