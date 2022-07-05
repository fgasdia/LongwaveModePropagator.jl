function test_physicalmodeequation(scenario)
    @unpack ea, tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    me = PhysicalModeEquation(ea, tx.frequency, waveguide)
    @test me isa LMP.ModeEquation
    @test me isa PhysicalModeEquation

    me2 = PhysicalModeEquation(tx.frequency, waveguide)
    @test me2 isa PhysicalModeEquation
    @test iszero(me2.θ)

    @test LMP.setea(ea, me2) == me
end

function test_wmatrix(scenario)
    @unpack ea, tx, bfield, species = scenario

    C = cos(ea)
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
            Wfcn(θ) = (T = LMP.tmatrix(θ, M); LMP.wmatrix(θ, T)[i][j])
            dWref = FiniteDiff.finite_difference_derivative(Wfcn, θs, Val{:central})
            dW(θ) = (T = LMP.tmatrix(θ, M); dT = LMP.dtmatrix(θ, M); LMP.dwmatrix(θ, T, dT)[i][j])

            @test maxabsdiff(dW.(θs), dWref) < 1e-3
        end
    end
end

function test_dRdz(scenario)
    @unpack ea, tx, bfield, species, ground = scenario

    params = LMPParams()
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)

    Mtop = LMP.susceptibility(params.topheight, me; params=params)
    Rtop = LMP.bookerreflection(ea, Mtop)

    # sharply bounded R from bookerreflection satisfies dR/dz = 0
    @test isapprox(LMP.dRdz(Rtop, me, params.topheight), zeros(2, 2); atol=1e-15)

    # Compare default and optional argument
    M = LMP.susceptibility(72e3, me; params=params)
    R = LMP.bookerreflection(ea, M)
    @inferred LMP.dRdz(R, me, 72e3)
end

function test_dRdθdz(scenario)
    @unpack ea, tx, bfield, species, ground = scenario

    params = LMPParams()
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)

    Mtop = LMP.susceptibility(params.topheight, me)
    z = params.topheight - 500

    for i = 1:4
        dRfcn(θ) = (Rtop = LMP.bookerreflection(θ, Mtop); me = LMP.setea(θ, me);
                LMP.dRdz(Rtop, me, z)[i])
        dRdθref = FiniteDiff.finite_difference_derivative(dRfcn, θs, Val{:central})
        dRdθtmp(θ) = (me = LMP.setea(θ, me);
                (Rtop, dRdθtop) = LMP.bookerreflection(θ, Mtop, LMP.Dθ());
                RdRdθtop = vcat(Rtop, dRdθtop);
                LMP.dRdθdz(RdRdθtop, (me, params), z))
        dR(θ) = dRdθtmp(θ)[SVector(1,2),:][i]
        dRdθ(θ) = dRdθtmp(θ)[SVector(3,4),:][i]

        @test dRfcn.(θs) ≈ dR.(θs)
        @test maxabsdiff(dRdθ.(θs), dRdθref) < 1e-6
    end
end

function test_integratedreflection_vertical(scenario)
    @unpack tx, ground, bfield, species = scenario

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(tx.frequency, waveguide)

    R = LMP.integratedreflection(me)

    @test R[1,2] ≈ R[2,1]
end

function test_integratedreflection_deriv(scenario)
    @unpack tx, ground, bfield, species = scenario
    freq = tx.frequency

    # 1e-10 is hardcoded in Dθ form of `integratedreflection`
    # This function is a test of the Dθ form - the lower default tolerance in the non-Dθ
    # form would cause tests to fail if not lowered to the same threshold
    ip = IntegrationParams(tolerance=1e-10)
    params = LMPParams(integrationparams=ip)
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(freq, waveguide)

    Rref(θ) = LMP.integratedreflection(LMP.setea(θ, me); params=params)
    RdR(θ) = LMP.integratedreflection(LMP.setea(θ, me), LMP.Dθ(); params=params)

    Rs = Vector{SMatrix{2,2,ComplexF64,4}}(undef, length(θs))
    dRs = similar(Rs)
    Rrefs = similar(Rs)
    dRrefs = similar(Rs)
    
    Threads.@threads for i in eachindex(θs)
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
        @test isapprox(R, Rr; rtol=1e-5)
        @test isapprox(dR, dRr; rtol=1e-3)
    end
end

function test_fresnelreflection(scenario)
    @unpack ea, tx, bfield, species, ground = scenario

    # PEC ground
    pec_ground = LMP.Ground(1, 1e12)
    vertical_ea = π/2
    Rg = LMP.fresnelreflection(vertical_ea, pec_ground, 24e3)
    @test isapprox(abs.(Rg), I; atol=1e-7)

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)
    @test LMP.fresnelreflection(ea, ground, tx.frequency) == LMP.fresnelreflection(me)
end

function test_fresnelreflection_deriv(scenario)
    @unpack ea, tx, bfield, species, ground = scenario
    freq = tx.frequency

    dRgtmp(θ) = LMP.fresnelreflection(θ, ground, freq, LMP.Dθ())
    for i = 1:4
        Rgref(θ) = LMP.fresnelreflection(θ, ground, freq)[i]
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

function test_solvemodalequation(scenario)
    @unpack ea, tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)

    f = LMP.solvemodalequation(me)
    me2 = PhysicalModeEquation(0.0+0.0im, tx.frequency, waveguide)
    f2 = LMP.solvemodalequation(ea, me2)
    @test f == f2
end

function test_modalequation_resonant(scenario)
    @unpack ea, tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    me = PhysicalModeEquation(ea, tx.frequency, waveguide)

    R = @SMatrix [1 0; 0 1]
    Rg = @SMatrix [1 0; 0 -1]
    @test isapprox(LMP.modalequation(R, Rg), 0; atol=1e-15)

    R = LMP.integratedreflection(me)
    Rg = LMP.fresnelreflection(ea, ground, tx.frequency)
    f = LMP.modalequation(R, Rg)
    @test isroot(f; atol=1e-4)  # this test is a little cyclic

    @test f == LMP.solvemodalequation(me)

    θ = 1.5 - 0.02im  # not resonant
    @test abs(LMP.solvemodalequation(θ, me)) > 1e-3
end

function test_modalequation_deriv(scenario)
    @unpack ea, tx, ground, bfield, species = scenario
    freq = tx.frequency

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    modeequation = PhysicalModeEquation(ea, freq, waveguide)

    # `solvedmodalequation` calls the `Dθ` form of `integratedreflection`, which uses a
    # hard coded tolerance of 1e-10. Tests will fail if comparing to the lower default
    # tolerance used by `solvemodalequation`
    ip = IntegrationParams(tolerance=1e-10)
    params = LMPParams(integrationparams=ip)

    dFref = FiniteDiff.finite_difference_derivative(θ->LMP.solvemodalequation(θ, modeequation; params=params),
        θs, Val{:central})

    dFs = Vector{ComplexF64}(undef, length(θs))
    Threads.@threads for i in eachindex(θs)
        dFdθ, R, Rg = LMP.solvedmodalequation(θs[i], modeequation)
        dFs[i] = dFdθ
    end

    @test isapprox(dFs, dFref; rtol=1e-3)
end

function test_findmodes(scenario)
    @unpack tx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    modeequation = PhysicalModeEquation(tx.frequency, waveguide)

    origcoords = LMP.defaultmesh(tx.frequency)

    modes = @inferred findmodes(modeequation, origcoords; params=LMPParams(refineeigenangles=false))
    modes2 = findmodes(modeequation, origcoords; params=LMPParams(refineeigenangles=true))
    modes3 = @inferred findmodes(modeequation)
    modes4 = findmodes(modeequation)

    @test modes2 ≈ modes4
    @test modes3 ≈ modes4
    @test length(modes) == length(modes2) == length(modes3)

    fparams = LMPParams(integrationparams=IntegrationParams(tolerance=1e-8))
    for m in modes
        f = LMP.solvemodalequation(m, modeequation; params=fparams)
        isroot(f; atol=1e-2) || return f
        @test isroot(f; atol=1e-2)
    end
    # fs = LMP.solvemodalequation.(modes, (modeequation,); params=fparams)

    for m in modes2
        f = LMP.solvemodalequation(m, modeequation; params=fparams)
        isroot(f; atol=1e-4) || return f
        @test isroot(f; atol=1e-4)
    end
    # fs2 = LMP.solvemodalequation.(modes2, (modeequation,); params=fparams)
end

########
# `TEST_MODES` must be filled for tests below!

function test_roots(scenario)
    x = 0.0005
    z = complex(x, -x)

    # isroot Real and Complex
    @test isroot(x)
    @test isroot(x; atol=1e-6) == false
    @test isroot(z)
    @test isroot(z; atol=1e-6) == false

    z2 = complex(1, 0)
    z3 = complex(0, 1)
    @test isroot(z2) == false
    @test isroot(z3) == false
end


@testset "modefinder.jl" begin
    @info "Testing modefinder"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        multiplespecies_scenario)
        
        test_physicalmodeequation(scn)  # just test with one scenario?

        test_wmatrix(scn)
        test_wmatrix_deriv(scn)

        test_dRdz(scn)
        test_dRdθdz(scn)

        test_integratedreflection_deriv(scn)

        test_fresnelreflection(scn)
        test_fresnelreflection_deriv(scn)

        test_solvemodalequation(scn)
        test_modalequation_deriv(scn)

        test_findmodes(scn)
    end

    # Fill in TEST_MODES
    @info "  Mode finding..."
    for scenario in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        multiplespecies_scenario)
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
