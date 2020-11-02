function randomwavefields()
    modes = EigenAngle.(rand(ComplexF64, 10))

    wavefields = LWMS.Wavefields(modes)

    return wavefields
end

function wavefields_test(scenario)
    @unpack tx, bfield, species, ground = scenario
    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide)

    # Δr from 0.5->0.25 => time from 3.8->5.3 sec
    # tolerance from 1e-8->1e-7 => time from 5.3->4.6 sec
    origcoords = LWMS.defaultcoordinates(tx.frequency)
    modes = LWMS.findmodes(modeequation, origcoords)

    wavefields = LWMS.Wavefields(modes)

    return isvalid(wavefields)
end

"""
Check for equivalence of Booker quartic computed by M and T.

Runtime is dominated by `roots!`, but the version with `T` is slightly faster.
"""
function booker_MTequivalence_test(scenario)
    @unpack ea, bfield, tx, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    qM, BM = LWMS.bookerquartic(ea, M)
    qT, BT = LWMS.bookerquartic(T)

    return sort(qM, by=real) ≈ sort(qT, by=real)
end

"""
Confirm bookerquartic!(T) finds valid eigenvalues.
"""
function booker_Tvalidity_test(scenario)
    @unpack ea, bfield, tx, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    qT, BT = LWMS.bookerquartic(T)

    # eigvals is >20 times slower than bookerquartic
    sort(eigvals(Array(T)), by=real) ≈ sort(qT, by=real) || return false

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(qT)
        booker = qT[i]^4 + BT[4]*qT[i]^3 + BT[3]*qT[i]^2 + BT[2]*qT[i] + BT[1]
        LWMS.isroot(booker, atol=1e-7) || return false
    end

    return true
end

"""
Confirm `initialwavefields(T)` satisfies field eigenvector equation ``Te = qe``
"""
function initialwavefields_test(scenario)
    @unpack ea, bfield, tx, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
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

    Mtop = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    Ttop = LWMS.tmatrix(ea, Mtop)
    etop = LWMS.initialwavefields(Ttop)

    wavefieldR = LWMS.vacuumreflectioncoeffs(ea, etop[:,1], etop[:,2])

    LWMS.sharpboundaryreflection(ea, Mtop) ≈ wavefieldR || return false

    return true
end


"""
Confirm reflection coefficients from wavefields match with dr/dz calculation.
"""
function drdzwavefield_equivalence_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario

    zs = LWMS.WAVEFIELD_HEIGHTS

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, species)
    wavefieldRs = [LWMS.vacuumreflectioncoeffs(ea, s[:,1], s[:,2]) for s in e]

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(ea, tx.frequency, waveguide)

    @unpack tolerance, solver = LWMS.DEFAULT_INTEGRATIONPARAMS

    Mtop = LWMS.susceptibility(first(zs), tx.frequency, waveguide)
    Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
    prob = ODEProblem{false}(LWMS.dRdz, Rtop, (first(zs), last(zs)), modeequation)
    sol = solve(prob, solver, abstol=tolerance, reltol=tolerance, # lower tolerance doesn't help match
                saveat=zs, save_everystep=false)

    return all(isapprox.(wavefieldRs, sol.u, rtol=1e-2))
end

"""
Check wavefields in homogeneous ionosphere are valid solutions to wave equation.

Compares to Booker quartic solution.

See, e.g. Pitteway 1965 pg 234; also Barron & Budden 1959 sec 10
"""
function homogeneous_iono_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario

    ionobottom = LWMS.CURVATURE_HEIGHT
    zs = 200e3:-500:ionobottom

    LWMS.set_earthcurvature(false)

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, species)

    # Normalize fields so component 2 (Ey) = 1, as is used in Booker Quartic
    e1 = [s[:,1]/s[2,1] for s in e]
    e2 = [s[:,2]/s[2,2] for s in e]

    e1 = reshape(reinterpret(ComplexF64, e1), 4, :)
    e2 = reshape(reinterpret(ComplexF64, e2), 4, :)

    # Booker solution - single solution for entire homogeneous iono
    M = LWMS.susceptibility(ionobottom, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)
    booker = LWMS.initialwavefields(T)

    all(e1 .≈ booker[:,1]) || return false
    all(e2 .≈ booker[:,2]) || return false

    # This is basically the same test...
    q, B = LWMS.bookerquartic(T)
    LWMS.sortquarticroots!(q)

    T*e1 ≈ q[1]*e1 || return false
    T*e2 ≈ q[2]*e2 || return false

    LWMS.set_earthcurvature(true)

    return true
end

function resonance_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario

    ztop = LWMS.TOPHEIGHT
    zs = ztop:-100:zero(LWMS.TOPHEIGHT)

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, species)
    R = LWMS.vacuumreflectioncoeffs(ea, e[end])
    Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)

    # Ensure we are close to mode resonance with R
    f = LWMS.modalequation(R, Rg)
    return LWMS.isroot(f)
end

function fieldstrengths_test(scenario)
    @unpack ea, bfield, tx, ground, species = scenario
    wavefields = randomwavefields()

    modes = LWMS.eigenangles(wavefields)

    LWMS.fieldstrengths!(view(wavefields,:,1), LWMS.WAVEFIELD_HEIGHTS, modes[1], tx.frequency,
                        bfield, species, ground)
end

@testset "wavefields.jl" begin
    @info "Testing wavefields"

    @test size(randomwavefields(),1) == length(LWMS.WAVEFIELD_HEIGHTS)
    @test length(LWMS.eigenangles(randomwavefields())) == 10

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        @test wavefields_test(scn)
        @test_nowarn fieldstrengths_test(scn)
    end

    @testset "Initial conditions" begin
        @info "  Initial conditions..."

        for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
            @test booker_MTequivalence_test(scn)
            @test booker_Tvalidity_test(scn)
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
            @test resonance_test(scn)
        end
        for scn in (homogeneousiono_scenario, )
            @test homogeneous_iono_test(scn)
        end
    end

    # TODO: Check that reflection coeffs for N/S directions are equal
    @test_skip wavefieldsR_north ≈ wavefieldsR_south

    # TODO: test `fieldstrengths`
    @test_skip LWMS.fieldstrengths()
end
