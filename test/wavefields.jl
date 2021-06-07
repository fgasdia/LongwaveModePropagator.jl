function randomwavefields()
    modes = EigenAngle.(rand(TEST_RNG, ComplexF64, 10))
    heights = LMPParams().wavefieldheights

    v = rand(SVector{6,ComplexF64}, length(heights), length(modes))
    wavefields = LMP.Wavefields(v, heights, modes)

    return modes, wavefields
end

function test_Wavefields()
    modes, wavefields = randomwavefields()
    @unpack wavefieldheights = LMPParams()

    @test LMP.numheights(wavefields) == length(wavefieldheights)
    @test LMP.nummodes(wavefields) == 10

    @test size(wavefields) == (length(wavefieldheights), 10)
    @test view(wavefields, 1:5) == view(wavefields.v, 1:5)

    newwavefields = similar(wavefields)
    @test size(newwavefields) == size(wavefields)
    @test LMP.heights(newwavefields) == LMP.heights(wavefields)
    @test LMP.eigenangles(newwavefields) == LMP.eigenangles(wavefields)

    newwavefields[1] = wavefields[1]*2
    @test newwavefields[1] == wavefields[1]*2

    # copy and ==
    cpwavefields = copy(wavefields)
    @test cpwavefields.heights == wavefields.heights
    @test cpwavefields.eas == wavefields.eas
    @test cpwavefields.v == wavefields.v
    @test cpwavefields == wavefields

    # isvalid
    heights = LMP.heights(wavefields)
    eas = LMP.eigenangles(wavefields)

    @test isvalid(wavefields)
    longerv = vcat(wavefields.v, rand(TEST_RNG, SVector{6,ComplexF64}, 1, length(eas)))
    badwavefields = LMP.Wavefields(longerv, heights, eas)
    @test !isvalid(badwavefields)

    widerv = hcat(wavefields.v, rand(TEST_RNG, SVector{6,ComplexF64}, length(heights)))
    badwavefields = LMP.Wavefields(widerv, heights, eas)
    @test !isvalid(badwavefields)
end

function test_WavefieldIntegrationParams(scenario)
    @unpack ea, tx, bfield, species = scenario

    params = LMPParams()
    topheight = params.topheight

    w1 = LMP.WavefieldIntegrationParams(topheight, LMP.BOTTOMHEIGHT, 0.0+0.0im, 1.0, 1.0, ea,
        tx.frequency, bfield, species, params)
    w2 = LMP.WavefieldIntegrationParams(topheight, ea, tx.frequency, bfield, species, params)

    @test w1 == w2
end

function test_integratewavefields_homogeneous(scenario)
    #==
    Check wavefields in homogeneous ionosphere are valid solutions to wave equation.
    Compares to Booker quartic solution.
    See, e.g. Pitteway 1965 pg 234; also Barron & Budden 1959 sec 10
    ==#

    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams(earthcurvature=false)

    ionobottom = params.curvatureheight
    zs = 200e3:-500:ionobottom

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species; params=params)

    # Normalize fields so component 2 (Ey) = 1, as is used in Booker Quartic
    e1 = [s[:,1]/s[2,1] for s in e]
    e2 = [s[:,2]/s[2,2] for s in e]

    e1 = reshape(reinterpret(ComplexF64, e1), 4, :)
    e2 = reshape(reinterpret(ComplexF64, e2), 4, :)

    # Booker solution - single solution for entire homogeneous iono
    M = LMP.susceptibility(ionobottom, tx.frequency, bfield, species; params=params)
    booker = LMP.bookerwavefields(ea, M)

    @test isapprox(e1[:,end], booker[:,1]; atol=1e-6)
    @test isapprox(e2[:,end], booker[:,2]; atol=1e-6)
end

function test_wavefieldreflection(scenario)
    # Confirm reflection coefficients from wavefields match with dr/dz calculation.

    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams()
    waveguide = HomogeneousWaveguide(bfield, species, ground)
    modeequation = PhysicalModeEquation(ea, tx.frequency, waveguide)

    zs = params.wavefieldheights

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species; params=params)
    wavefieldRs = [LMP.bookerreflection(ea, s) for s in e]

    @unpack tolerance, solver = params.integrationparams

    Mtop = LMP.susceptibility(first(zs), tx.frequency, waveguide)
    Rtop = LMP.bookerreflection(ea, Mtop)
    prob = ODEProblem{false}(LMP.dRdz, Rtop, (first(zs), last(zs)), modeequation)
    sol = solve(prob, solver; abstol=tolerance, reltol=tolerance,
                saveat=zs, save_everystep=false)

    @test isapprox(wavefieldRs, sol.u, rtol=1e-3)
end

function test_wavefieldreflection_resonant(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams()

    ztop = params.topheight
    zs = ztop:-100:0.0

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species; params=params)
    R = LMP.bookerreflection(ea, e[end])
    Rg = LMP.fresnelreflection(ea, ground, tx.frequency)

    # Ensure we are close to mode resonance with R
    f = LMP.modalequation(R, Rg)

    @test isroot(f)
end

function test_boundaryscalars(scenario)
    # Check if boundaryscalars match for isotropic=false and isotropic=true when both
    # actually are isotropic

    @unpack ea, bfield, tx, ground, species = scenario
    params = LMPParams()

    ztop = params.topheight
    zs = ztop:-100:0.0

    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species; params=params)
    R = LMP.bookerreflection(ea, e[end])
    Rg = LMP.fresnelreflection(ea, ground, tx.frequency)

    a1, a2 = LMP.boundaryscalars(R, Rg, e[end], false)
    b1, b2 = LMP.boundaryscalars(R, Rg, e[end], true)

    # return a1, a2, b1, b2

    @test a1 ≈ b1
    @test a2 ≈ b2
end

function test_fieldstrengths(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    modes, wavefields = randomwavefields()

    modes = LMP.eigenangles(wavefields)

    LMP.fieldstrengths!(view(wavefields,:,1), LMPParams().wavefieldheights, modes[1],
                         tx.frequency, bfield, species, ground)
end

@testset "wavefields.jl" begin
    @info "Testing wavefields"

    test_Wavefields()

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        test_WavefieldIntegrationParams(scn)
        @test_nowarn test_fieldstrengths(scn)
    end

    @testset "Integration" begin
        for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
            test_wavefieldreflection(scn)
        end
        for scn in (resonant_scenario, )
            test_wavefieldreflection_resonant(scn)
        end
        for scn in (homogeneousiono_scenario, )
            test_integratewavefields_homogeneous(scn)
        end

        @test_skip test_boundaryscalars(isotropicB_resonant_scenario)
    end
end
