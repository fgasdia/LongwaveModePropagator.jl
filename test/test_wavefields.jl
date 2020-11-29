function randomwavefields()
    modes = EigenAngle.(rand(TEST_RNG, ComplexF64, 10))

    wavefields = LWMS.Wavefields(LWMSParams().wavefieldheights, modes)

    return wavefields
end

function wavefieldfuncs_test()
    modes = EigenAngle.(rand(TEST_RNG, ComplexF64, 10))

    @unpack wavefieldheights = LWMSParams()

    wavefields = LWMS.Wavefields(wavefieldheights, modes)

    LWMS.numheights(wavefields) == length(wavefieldheights) || return false
    LWMS.nummodes(wavefields) == 10 || return false

    return true
end

function validwavefields_test()
    wavefields = randomwavefields()

    heights = LWMS.heights(wavefields)
    eas = LWMS.eigenangles(wavefields)

    isvalid(wavefields) || return false

    longerv = vcat(wavefields.v, rand(TEST_RNG, SVector{6,ComplexF64}, 1, length(eas)))
    wavefields = LWMS.Wavefields(longerv, heights, eas)
    !isvalid(wavefields) || return false

    wavefields = randomwavefields()
    widerv = hcat(wavefields.v, rand(TEST_RNG, SVector{6,ComplexF64}, length(heights)))
    wavefields = LWMS.Wavefields(widerv, heights, eas)
    !isvalid(wavefields) || return false

    return true
end

"""
Confirm `initialwavefields(T)` satisfies field eigenvector equation ``Te = qe``
"""
function initialwavefields_test(scenario)
    @unpack ea, bfield, tx, species = scenario

    topheight = first(LWMSParams().wavefieldheights)
    M = LWMS.susceptibility(topheight, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    # Verify initialwavefields produces a valid solution
    e = LWMS.initialwavefields(T)
    q, B = LWMS.bookerquartic(T)
    LWMS.sortquarticroots!(q)

    for i = 1:2
        T*e[:,i] ≈ q[i]*e[:,i] || return false
    end

    return true
end

"""
Confirm `vacuumreflectioncoeffs` agrees with `sharpboundaryreflection` at the
top height.
"""
function initialR_test(scenario)
    @unpack ea, bfield, tx, species = scenario
    params = LWMSParams()

    topheight = first(params.wavefieldheights)
    Mtop = LWMS.susceptibility(topheight, tx.frequency, bfield, species, params=params)
    Ttop = LWMS.tmatrix(ea, Mtop)
    etop = LWMS.initialwavefields(Ttop)

    wavefieldR = LWMS.vacuumreflectioncoeffs(ea, etop[:,1], etop[:,2])

    Rref = LWMS.sharpboundaryreflection(ea, Mtop)
    isapprox(Rref, wavefieldR, rtol=1e-6) || return false

    return true
end


"""
Confirm reflection coefficients from wavefields match with dr/dz calculation.
"""
function drdzwavefield_equivalence_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    params = LWMSParams()

    zs = params.wavefieldheights

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, species)

    wavefieldRs = [LWMS.vacuumreflectioncoeffs(ea, s[:,1], s[:,2]) for s in e]

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(ea, tx.frequency, waveguide)

    @unpack tolerance, solver = params.integrationparams

    Mtop = LWMS.susceptibility(first(zs), tx.frequency, waveguide)
    Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
    prob = ODEProblem{false}(LWMS.dRdz, Rtop, (first(zs), last(zs)), (modeequation, params))
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
    params = LWMSParams(earthcurvature=false)

    ionobottom = params.curvatureheight
    zs = 200e3:-500:ionobottom

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, species, params=params)

    # Normalize fields so component 2 (Ey) = 1, as is used in Booker Quartic
    e1 = [s[:,1]/s[2,1] for s in e]
    e2 = [s[:,2]/s[2,2] for s in e]

    e1 = reshape(reinterpret(ComplexF64, e1), 4, :)
    e2 = reshape(reinterpret(ComplexF64, e2), 4, :)

    # Booker solution - single solution for entire homogeneous iono
    M = LWMS.susceptibility(ionobottom, tx.frequency, bfield, species, params=params)
    T = LWMS.tmatrix(ea, M)
    booker = LWMS.initialwavefields(T)

    isapprox(e1[:,end], booker[:,1], rtol=1e-6) || return false
    isapprox(e2[:,end], booker[:,2], rtol=1e-6) || return false

    # This is basically the same test...
    q, B = LWMS.bookerquartic(T)
    LWMS.sortquarticroots!(q)

    isapprox(T*e1[:,end], q[1]*e1[:,end], rtol=1e-6) || return false
    isapprox(T*e2[:,end], q[2]*e2[:,end], rtol=1e-6) || return false

    return true
end

function resonance_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    params = LWMSParams()https://github.com/fgasdia/RootsAndPoles.jl

    ztop = params.topheight
    zs = ztop:-100:0.0


    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide)

    origcoords = LWMS.defaultcoordinates(tx.frequency)
    modes = LWMS.findmodes(modeequation, origcoords, params=params)
    return modes
    ea = modes[1]

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, species, params=params)
    R = LWMS.vacuumreflectioncoeffs(ea, e[end])
    Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)

    # Ensure we are close to mode resonance with R
    f = LWMS.modalequation(R, Rg)

    return LWMS.isroot(f)
end

function boundary_test(scenario)
    # TODO: Create a vertical resonant scenario
    # Check if boundaryscalars match for isotropic=false and isotropic=true when both
    # actually are isotropic

    @unpack ea, bfield, tx, ground, species = scenario
    params = LWMSParams()

    ztop = params.topheight
    zs = ztop:-100:0.0

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, species, params=params)
    R = LWMS.vacuumreflectioncoeffs(ea, e[end])
    Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)

    a1, a2 = LWMS.boundaryscalars(R, Rg, e[end], false)
    b1, b2 = LWMS.boundaryscalars(R, Rg, e[end], true)

    # return a1, a2, b1, b2

    a1 ≈ b1 || return false
    a2 ≈ b2 || return false

    return true
end

function fieldstrengths_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    wavefields = randomwavefields()

    modes = LWMS.eigenangles(wavefields)

    LWMS.fieldstrengths!(view(wavefields,:,1), LWMSParams().wavefieldheights, modes[1],
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
            @test initialwavefields_test(scn)
            @test initialR_test(scn)
        end
    end

    @testset "Integration" begin
        @info "  Wavefield integration..."

        for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
            @test drdzwavefield_equivalence_test(scn)
        end
        for scn in (resonant_scenario, )
            @test_broken resonance_test(scn)
        end
        for scn in (homogeneousiono_scenario, )
            @test homogeneous_iono_test(scn)
        end

        @test_broken boundary_test(isotropicB_resonant_scenario)
    end

    # TODO: test `fieldstrengths`, but how?
end
