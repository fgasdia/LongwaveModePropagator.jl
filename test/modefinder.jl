function test_physicalmodeequation(scenario)
    @unpack ea, tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    me = PhysicalModeEquation(ea, tx.frequency, waveguide)
    @test me isa LMP.ModeEquation
    @test me isa PhysicalModeEquation

    me2 = PhysicalModeEquation(tx.frequency, waveguide)
    @test me2 isa PhysicalModeEquation
    @test iszero(me2.ea.θ)

    @test LMP.setea(ea, me2) == me
    @test LMP.setea(ea.θ, me2) == me
end

function test_roots(scenario)
    x = 0.0005
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
    @unpack ea, tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(tx.frequency, waveguide)

    # nothing filtered
    roots = copy(TEST_MODES[scenario])
    @test LMP.filterroots!(roots, me) == roots
    @test LMP.filterroots!(roots, tx.frequency, waveguide) == roots

    # filter bad root
    push!(roots, EigenAngle(complex(1.5, -0.5)))
    @test LMP.filterroots!(roots, me) == TEST_MODES[scenario]
    @test LMP.filterroots!(roots, tx.frequency, waveguide) == TEST_MODES[scenario]

    # filter (or not) with different tolerance
    roots2 = copy(roots)
    roots2[1] = EigenAngle(roots2[1].θ + 0.001)
    f = LMP.solvemodalequation(LMP.setea(roots2[1], me))
    @test LMP.isroot(f, atol=0.5)  # NOTE: atol is for value of modal equation, not θ
    @test LMP.filterroots!(roots2, me, atol=0.5) == roots2
    @test LMP.filterroots!(roots2, me) == roots[2:end]
end

function test_wmatrix(scenario)
    @unpack ea, tx, bfield, species = scenario

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
    @unpack tx, bfield, species = scenario

    M = LMP.susceptibility(80e3, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Wfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M); LMP.wmatrix(ea, T)[i][j])
            dWref = FiniteDiff.finite_difference_derivative(Wfcn, θs, Val{:central})
            dW(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M);
                     dT = LMP.dtmatrix(ea, M); LMP.dwmatrix(ea, T, dT)[i][j])

            @test maxabsdiff(dW.(θs), dWref) < 1e-3
        end
    end
end

function test_dRdz(scenario)
    @unpack ea, tx, bfield, species, ground = scenario

    params = LMPParams()
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)

    Mtop = LMP.susceptibility(params.topheight, me, params=params)
    Rtop = LMP.bookerreflection(ea, Mtop)

    # sharply bounded R from bookerreflection satisfies dR/dz = 0
    @test isapprox(LMP.dRdz(Rtop, (me, params), params.topheight), zeros(2, 2), atol=1e-15)
end

function test_dRdθdz(scenario)
    @unpack ea, tx, bfield, species, ground = scenario

    params = LMPParams()
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)

    Mtop = LMP.susceptibility(params.topheight, me)
    z = params.topheight - 500

    for i = 1:4
        dRfcn(θ) = (ea = EigenAngle(θ); Rtop = LMP.bookerreflection(ea, Mtop);
            me = LMP.setea(ea, me); LMP.dRdz(Rtop, (me, params), z)[i])
        dRdθref = FiniteDiff.finite_difference_derivative(dRfcn, θs, Val{:central})
        dRdθtmp(θ) = (ea = EigenAngle(θ); me = LMP.setea(ea, me);
            (Rtop, dRdθtop) = LMP.bookerreflection(ea, Mtop, LMP.Dθ());
            RdRdθtop = vcat(Rtop, dRdθtop);
            LMP.dRdθdz(RdRdθtop, (me, params), z))
        dR(θ) = dRdθtmp(θ)[SVector(1,2),:][i]
        dRdθ(θ) = dRdθtmp(θ)[SVector(3,4),:][i]

        @test dRfcn.(θs) ≈ dR.(θs)
        @test maxabsdiff(dRdθ.(θs), dRdθref) < 1e-6
    end
end

function test_integratedreflection_vertical(scenario)
    @unpack ea, tx, ground, bfield, species = scenario

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(tx.frequency, waveguide)

    R = LMP.integratedreflection(me)

    @test R[1,2] ≈ R[2,1]
end

function test_integratedreflection_deriv(scenario)
    @unpack tx, ground, bfield, species = scenario
    freq = tx.frequency

    params = LMPParams()
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(freq, waveguide)

    Rref(θ) = (me = LMP.setea(θ, me); LMP.integratedreflection(me, params=params))
    RdR(θ) = (me = LMP.setea(θ, me); LMP.integratedreflection(me, LMP.Dθ(), params=params))

    Rs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(θs))
    dRs = similar(Rs)
    Rrefs = similar(Rs)
    dRrefs = similar(Rs)
    # NOTE: this also effectively checks for thread safety
    Threads.@threads for i in 1:length(θs)
        v = RdR(θs[i])
        Rs[i] = v[SVector(1,2),:]
        dRs[i] = v[SVector(3,4),:]
        Rrefs[i] = Rref(θs[i])
        dRrefs[i] = FiniteDiff.finite_difference_derivative(Rref, θs[i], Val{:central})
    end

    for i = 1:4
        R = getindex.(Rs, i)
        dR = getindex.(dRs, i)
        Rr = getindex.(Rrefs, i)
        dRr = getindex.(dRrefs, i)

        # maxabsdiff criteria doesn't capture range of R and dR so rtol is used
        @test isapprox(R, Rr, rtol=1e-5)
        @test isapprox(dR, dRr, rtol=1e-3)
    end
end

function test_fresnelreflection(scenario)
    @unpack ea, tx, bfield, species, ground = scenario

    # PEC ground
    pec_ground = LMP.Ground(1, 1e12)
    vertical_ea = LMP.EigenAngle(π/2)
    Rg = LMP.fresnelreflection(vertical_ea, pec_ground, Frequency(24e3))
    @test isapprox(abs.(Rg), I, atol=1e-7)

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)
    @test LMP.fresnelreflection(ea, ground, tx.frequency) == LMP.fresnelreflection(me)
end

function test_fresnelreflection_deriv(scenario)
    @unpack ea, tx, bfield, species, ground = scenario
    freq = tx.frequency

    dRgtmp(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq, LMP.Dθ()))
    for i = 1:4
        Rgref(θ) = (ea = EigenAngle(θ); LMP.fresnelreflection(ea, ground, freq)[i])
        dRgref = FiniteDiff.finite_difference_derivative(Rgref, θs, Val{:central})
        Rg(θ) = dRgtmp(θ)[1][i]
        dRg(θ) = dRgtmp(θ)[2][i]

        @test Rgref.(θs) ≈ Rg.(θs)
        @test maxabsdiff(dRg.(θs), dRgref) < 1e-6
    end

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)
    Rg1 = LMP.fresnelreflection(ea, ground, freq, LMP.Dθ())
    Rg2 = LMP.fresnelreflection(me, LMP.Dθ())
    @test Rg1 == Rg2
end

function test_modalequation_resonant(scenario)
    @unpack ea, tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)

    R = @SMatrix [1 0; 0 1]
    Rg = @SMatrix [1 0; 0 -1]
    @test isapprox(LMP.modalequation(R, Rg), 0, atol=1e-15)

    R = LMP.integratedreflection(me)
    Rg = LMP.fresnelreflection(ea, ground, tx.frequency)
    f = LMP.modalequation(R, Rg)
    @test LMP.isroot(f, atol=1e-4)  # this test is a little cyclic

    @test f == LMP.solvemodalequation(me)

    θ = 1.5 - 0.02im  # not resonant
    @test abs(LMP.solvemodalequation(θ, me)) > 1e-3
end

function test_modalequation_deriv(scenario)
    @unpack ea, tx, ground, bfield, species = scenario
    freq = tx.frequency

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    modeequation = PhysicalModeEquation(ea, freq, waveguide)

    dFref = FiniteDiff.finite_difference_derivative(θ->LMP.solvemodalequation(θ, modeequation),
        θs, Val{:central})

    dFs = Vector{ComplexF64}(undef, length(θs))
    Threads.@threads for i in eachindex(θs)
        dFdθ, R, Rg = LMP.solvedmodalequation(θs[i], modeequation)
        dFs[i] = dFdθ
    end

    @test isapprox(dFs, dFref, rtol=1e-3)
end

function test_findmodes(scenario)
    @unpack tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    modeequation = PhysicalModeEquation(tx.frequency, waveguide)

    origcoords = LMP.defaultmesh(tx.frequency)

    # params = LMPParams(grpfparams=LMP.GRPFParams(100000, 1e-6, true))
    params = LMPParams()
    modes = findmodes(modeequation, origcoords, params=params)

    for m in modes
        f = LMP.solvemodalequation(m, modeequation, params=params)
        LMP.isroot(f) || return f
        @test LMP.isroot(f)
    end

    # return modes
end

# function evalroot(root, scenario)
#     @unpack tx, bfield, species, ground = scenario
#     waveguide = HomogeneousWaveguide(bfield, species, ground)
#
#     modeequation = PhysicalModeEquation(tx.frequency, waveguide)
#     LMP.solvemodalequation(EigenAngle(root), modeequation)
# end


@testset "modefinder.jl" begin
    @info "Testing modefinder"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        test_physicalmodeequation(scn)  # just test with one scenario?

        test_wmatrix(scn)
        test_wmatrix_deriv(scn)

        test_dRdz(scn)
        test_dRdθdz(scn)

        test_integratedreflection_deriv(scn)

        test_fresnelreflection(scn)
        test_fresnelreflection_deriv(scn)

        test_modalequation_deriv(scn)

        test_findmodes(scn)
    end

    # Fill in TEST_MODES
    @info "  Mode finding..."
    for scenario in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        if !haskey(TEST_MODES, scenario)
            TEST_MODES[scenario] = findroots(scenario)
        end
    end

    for scn in (resonant_scenario, )
        test_roots(scn)
        test_modalequation_resonant(scn)
    end
    for scn in (verticalB_scenario, )
        test_integratedreflection_vertical(scn)
    end
end
