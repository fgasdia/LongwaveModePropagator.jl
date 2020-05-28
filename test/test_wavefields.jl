function scenario()
    bfield = BField(50e-6, deg2rad(68), deg2rad(111))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
    ground = Ground(15, 0.001)

    electrons = Constituent(qₑ, mₑ,
                            z -> waitprofile(z, 75, 0.32),
                            electroncollisionfrequency)

    ea = EigenAngle(1.45964665843992 - 0.014974434753336im)

    ztop = LWMS.TOPHEIGHT
    zs = ztop:-100:zero(LWMS.TOPHEIGHT)

    return bfield, tx, ground, electrons, ea, zs
end


"""
Check for equivalence of Booker quartic computed by M and T.

Runtime is dominated by `roots!`, but the version with `T` is slightly faster.
"""
function booker_MTequivalence_test()
    bfield, tx, ground, electrons, ea, zs = scenario()

    M = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)

    LWMS.bookerquartic!(ea, M)
    qM, BM = copy(LWMS.BOOKER_QUARTIC_ROOTS), copy(LWMS.BOOKER_QUARTIC_COEFFS)

    LWMS.bookerquartic!(T)
    qT, BT = copy(LWMS.BOOKER_QUARTIC_ROOTS), copy(LWMS.BOOKER_QUARTIC_COEFFS)

    return sort(qM, by=real) ≈ sort(qT, by=real)
end

"""
Confirm bookerquartic!(T) finds valid eigenvalues.
"""
function booker_Tvalidity_test()
    bfield, tx, ground, electrons, ea, zs = scenario()

    M = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)

    LWMS.bookerquartic!(T)
    qT, BT = copy(LWMS.BOOKER_QUARTIC_ROOTS), copy(LWMS.BOOKER_QUARTIC_COEFFS)

    # eigvals is >20 times slower than bookerquartic
    sort(eigvals(Array(T)), by=real) ≈ sort(qT, by=real) || return false

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(qT)
        booker = qT[i]^4 + BT[4]*qT[i]^3 + BT[3]*qT[i]^2 + BT[2]*qT[i] + BT[1]
        isapprox(booker, 0, atol=1e-7) || return false
    end

    return true
end

"""
Confirm `initialwavefields(T)` satisfies field eigenvector equation ``Te = qe``
"""
function initialwavefields_test()
    bfield, tx, ground, electrons, ea, zs = scenario()

    M = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)

    # Verify initialwavefields produces a valid solution
    e = LWMS.initialwavefields(T)
    q = LWMS.BOOKER_QUARTIC_ROOTS

    for i = 1:2
        T*e[:,i] ≈ q[i]*e[:,i] || return false
    end

    return true
end

"""
Confirm `vacuumreflectioncoeffs` agrees with `sharpboundaryreflection` at the
top height.
"""
function initialR_test()
    # Confirm reflection coefficient from wavefields at top height
    bfield, tx, ground, electrons, ea, zs = scenario()

    Mtop = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
    Ttop = LWMS.tmatrix(ea, Mtop)
    etop = LWMS.initialwavefields(Ttop)

    wavefieldR = LWMS.vacuumreflectioncoeffs(ea, etop[:,1], etop[:,2])

    # "modulus" (abs) of each component should be <=1
    all(abs.(wavefieldR) .<= 1) || return false
    LWMS.sharpboundaryreflection(ea, Mtop) ≈ wavefieldR || return false

    return true
end


"""
Confirm reflection coefficients from wavefields match with dr/dz calculation.
"""
function drdzwavefield_equivalence_test()
    bfield, tx, ground, electrons, ea, zs = scenario()

    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, electrons)
    wavefieldRs = [LWMS.vacuumreflectioncoeffs(ea, s[:,1], s[:,2]) for s in e]

    modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)
    Mtop = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
    Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
    prob = ODEProblem{false}(LWMS.dRdz, Rtop, (first(zs), last(zs)), (ea, modeparams))
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8,
                saveat=zs, save_everystep=false)

    return all(isapprox.(wavefieldRs, sol.u, atol=1e-7))
end

function homogeneous_scenario()
    bfield = BField(50e-6, deg2rad(68), deg2rad(111))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
    ground = Ground(15, 0.001)

    ionobottom = 50e3
    species = Constituent(qₑ, mₑ, z -> z >= ionobottom ? ionobottom : 0.0,
                              z -> 5e6)  # ν is constant

    # Resonant EigenAngle
    ea = EigenAngle(1.45964665843992 - 0.014974434753336im)

    zs = 200e3:-500:ionobottom

    return bfield, tx, ground, species, ea, zs
end

"""
Check wavefields in homogeneous ionosphere are valid solutions to wave equation.

Compares to Booker quartic solution.

See, e.g. Pitteway 1965 pg 234; also Barron & Budden 1959 sec 10
"""
function homogeneous_iono_test()
    bfield, tx, ground, species, ea, zs = homogeneous_scenario()

    LWMS.EARTHCURVATURE[] = false

    ionobottom = last(zs)

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

    e1diff = e1 .- booker[:,1]
    e2diff = e2 .- booker[:,2]

    q = LWMS.BOOKER_QUARTIC_ROOTS

    all(e1 .≈ booker[:,1]) || return false
    all(e2 .≈ booker[:,2]) || return false

    # This is basically the same test...
    T*e1 ≈ q[1]*e1 || return false
    T*e2 ≈ q[2]*e2 || return false

    LWMS.EARTHCURVATURE[] = true

    return true
end

function resonant_scenario()
    bfield, tx, ground, electrons, ea, zs = scenario()
    modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8

    modes = LWMS.findmodes(origcoords, modeparams, tolerance)
    ea = modes[argmax(real(modes))]  # largest real resonant mode

    return bfield, tx, ground, electrons, EigenAngle(ea), zs
end

function resonance_test()
    bfield, tx, ground, electrons, ea, zs = resonant_scenario()
    e = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, electrons)
    R = LWMS.vacuumreflectioncoeffs(ea, e[end])
    Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)
    b1, b2 = LWMS.boundaryscalars(R, Rg, e[end])

    # Ensure we are close to mode resonance with R
    f = LWMS.modalequation(R, Rg)
    isapprox(f, 0, atol=1e-6)
end

@testset "Wavefields" begin
    @info "Testing wavefield functions..."

    @testset "Initial conditions" begin
        @info "  Testing initial conditions..."

        @test booker_MTequivalence_test()
        @test booker_Tvalidity_test()
        @test initialwavefields_test()
        @test initialR_test()
    end

    @testset "Integration" begin
        @info "  Testing wavefield integration..."

        @test drdzwavefield_equivalence_test()
        @test homogeneous_iono_test()
        @test resonance_test()
    end

    # TODO: Check that reflection coeffs for N/S directions are equal
    @test_skip wavefieldsR_north ≈ wavefieldsR_south

    # TODO: test `fieldstrengths`
    @test_skip LWMS.fieldstrengths()
end
