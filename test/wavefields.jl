function randomwavefields()
    modes = EigenAngle.(rand(TEST_RNG, ComplexF64, 10))

    wavefields = LMP.Wavefields(LMPParams().wavefieldheights, modes)

    return wavefields
end

function wavefieldfuncs_test()
    modes = EigenAngle.(rand(TEST_RNG, ComplexF64, 10))

    @unpack wavefieldheights = LMPParams()

    wavefields = LMP.Wavefields(wavefieldheights, modes)

    LMP.numheights(wavefields) == length(wavefieldheights) || return false
    LMP.nummodes(wavefields) == 10 || return false

    return true
end

function validwavefields_test()
    wavefields = randomwavefields()

    heights = LMP.heights(wavefields)
    eas = LMP.eigenangles(wavefields)

    isvalid(wavefields) || return false

    longerv = vcat(wavefields.v, rand(TEST_RNG, SVector{6,ComplexF64}, 1, length(eas)))
    wavefields = LMP.Wavefields(longerv, heights, eas)
    !isvalid(wavefields) || return false

    wavefields = randomwavefields()
    widerv = hcat(wavefields.v, rand(TEST_RNG, SVector{6,ComplexF64}, length(heights)))
    wavefields = LMP.Wavefields(widerv, heights, eas)
    !isvalid(wavefields) || return false

    return true
end


"""
Confirm `vacuumreflectioncoeffs` agrees with `bookerreflection` at the
top height.
"""
function initialR_test(scenario)
    @unpack ea, bfield, tx, species = scenario
    params = LMPParams()

    topheight = first(params.wavefieldheights)
    Mtop = LMP.susceptibility(topheight, tx.frequency, bfield, species, params=params)
    Ttop = LMP.tmatrix(ea, Mtop)
    etop = LMP.bookerwavefields(Ttop)

    e1, e2 = etop[:,1], etop[:,2]

    wavefieldR = LMP.vacuumreflectioncoeffs(ea, e1, e2)

    Cinv = ea.secθ
    Sv_inv = SMatrix{4,4}(Cinv, 0, -Cinv, 0,
                          0, -1, 0, -1,
                          0, -Cinv, 0, Cinv,
                          1, 0, 1, 0)

    f1 = Sv_inv*e1
    f2 = Sv_inv*e2

    U = SMatrix{2,2}(f1[1], f1[2], f2[1], f2[2])
    D = SMatrix{2,2}(f1[3], f1[4], f2[3], f2[4])

    Rref = D/U
    Rref ≈ wavefieldR || return false

    Rref = LMP.bookerreflection(ea, Mtop)
    Rref ≈ wavefieldR || return false

    return true
end


"""
Confirm reflection coefficients from wavefields match with dr/dz calculation.
"""
function drdzwavefield_equivalence_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams()

    zs = params.wavefieldheights

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species)

    wavefieldRs = [LMP.vacuumreflectioncoeffs(ea, s[:,1], s[:,2]) for s in e]

    waveguide = LMP.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LMP.PhysicalModeEquation(ea, tx.frequency, waveguide)

    @unpack tolerance, solver = params.integrationparams

    Mtop = LMP.susceptibility(first(zs), tx.frequency, waveguide)
    Rtop = LMP.bookerreflection(ea, Mtop)
    prob = ODEProblem{false}(LMP.dRdz, Rtop, (first(zs), last(zs)), (modeequation, params))
    sol = solve(prob, solver, abstol=tolerance, reltol=tolerance, # lower tolerance doesn't help match
                saveat=zs, save_everystep=false)

    return all(isapprox.(wavefieldRs, sol.u, rtol=1e-3))
end

"""
Check wavefields in homogeneous ionosphere are valid solutions to wave equation.

Compares to Booker quartic solution.

See, e.g. Pitteway 1965 pg 234; also Barron & Budden 1959 sec 10
"""
function homogeneous_iono_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams(earthcurvature=false)

    ionobottom = params.curvatureheight
    zs = 200e3:-500:ionobottom

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species, params=params)

    # Normalize fields so component 2 (Ey) = 1, as is used in Booker Quartic
    e1 = [s[:,1]/s[2,1] for s in e]
    e2 = [s[:,2]/s[2,2] for s in e]

    e1 = reshape(reinterpret(ComplexF64, e1), 4, :)
    e2 = reshape(reinterpret(ComplexF64, e2), 4, :)

    # Booker solution - single solution for entire homogeneous iono
    M = LMP.susceptibility(ionobottom, tx.frequency, bfield, species, params=params)
    T = LMP.tmatrix(ea, M)
    booker = LMP.bookerwavefields(T)

    isapprox(e1[:,end], booker[:,1], rtol=1e-6) || return false
    isapprox(e2[:,end], booker[:,2], rtol=1e-6) || return false

    # This is basically the same test...
    q, B = LMP.bookerquartic(T)
    LMP.sortquarticroots!(q)

    isapprox(T*e1[:,end], q[1]*e1[:,end], rtol=1e-6) || return false
    isapprox(T*e2[:,end], q[2]*e2[:,end], rtol=1e-6) || return false

    return true
end

function resonance_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams()

    ztop = params.topheight
    zs = ztop:-100:0.0

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species, params=params)
    R = LMP.vacuumreflectioncoeffs(ea, e[end])
    Rg = LMP.fresnelreflection(ea, ground, tx.frequency)

    # Ensure we are close to mode resonance with R
    f = LMP.modalequation(R, Rg)

    return LMP.isroot(f)
end

function boundary_test(scenario)
    # TODO: Create a vertical resonant scenario
    # Check if boundaryscalars match for isotropic=false and isotropic=true when both
    # actually are isotropic

    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams()

    ztop = params.topheight
    zs = ztop:-100:0.0

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species, params=params)
    R = LMP.vacuumreflectioncoeffs(ea, e[end])
    Rg = LMP.fresnelreflection(ea, ground, tx.frequency)

    a1, a2 = LMP.boundaryscalars(R, Rg, e[end], false)
    b1, b2 = LMP.boundaryscalars(R, Rg, e[end], true)

    # return a1, a2, b1, b2

    a1 ≈ b1 || return false
    a2 ≈ b2 || return false

    return true
end

function fieldstrengths_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    wavefields = randomwavefields()

    modes = LMP.eigenangles(wavefields)

    LMP.fieldstrengths!(view(wavefields,:,1), LMPParams().wavefieldheights, modes[1],
                         tx.frequency, bfield, species, ground)
end

@testset "wavefields.jl" begin
    @info "Testing wavefields"

    @test wavefieldfuncs_test()
    @test validwavefields_test()

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        @test_nowarn fieldstrengths_test(scn)
    end

    @testset "Initial conditions" begin
        @info "  Initial conditions..."

        for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
            @test bookerwavefields_test(scn)
            @test initialR_test(scn)
        end
    end

    @testset "Integration" begin
        @info "  Wavefield integration..."

        for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
            @test drdzwavefield_equivalence_test(scn)
        end
        for scn in (resonant_scenario, )
            @test resonance_test(scn)
        end
        for scn in (homogeneousiono_scenario, )
            @test homogeneous_iono_test(scn)
        end

        @test_skip boundary_test(isotropicB_resonant_scenario)
    end

    # TODO: test `fieldstreng9ths`, but how?
end
