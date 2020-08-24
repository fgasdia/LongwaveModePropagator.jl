const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C

const Δθ = deg2rad(0.01)

centereddiff(v, Δθ) = (v[3:end] - v[1:end-2])/(2Δθ)  # aligns to [2:end-1]

function derivative_scenario()
    # NOTE: as a finite diff, Δθ determines the accuracy
    # Setting too small will _greatly_ increase test runtime
    eas = [EigenAngle(th) for th in complex.(range(0, π/2; step=Δθ), -4π/180)]

    freq = Frequency(24e3)
    bfield = BField(50e-6, π/2, 0)
    ground = Ground(15, 0.001)
    electrons = Species(qₑ, mₑ,
                        z -> waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
                        z -> electroncollisionfrequency(z, LWMS.CURVATURE_HEIGHT))

    return eas, freq, ground, bfield, electrons
end

function wmatrix_deriv()
    eas, freq, ground, bfield, electrons = derivative_scenario()

    M = LWMS.susceptibility(70e3, freq, bfield, electrons)

    Ws = NTuple{4,SArray{Tuple{2,2},ComplexF64,2,4}}[]
    dWs = NTuple{4,SArray{Tuple{2,2},ComplexF64,2,4}}[]
    for ea in eas
        T = LWMS.tmatrix(ea, M)
        W = LWMS.wmatrix(ea, T)
        dW = LWMS.dwmatrixdθ(ea, M, T)
        push!(Ws, W)
        push!(dWs, dW)
    end

    for i = 1:4
        # Extract each tuple matrix of Ws
        ws = [w[i] for w in Ws]
        dws = [dw[i] for dw in dWs]

        re_ws = Tuple.(real.(ws))
        im_ws = Tuple.(imag.(ws))

        re_dws = Tuple.(real.(dws))
        im_dws = Tuple.(imag.(dws))

        for j = 1:4
            re_diff = centereddiff(getfield.(re_ws, j), Δθ) - getfield.(re_dws, j)[2:end-1]
            im_diff = centereddiff(getfield.(im_ws, j), Δθ) - getfield.(im_dws, j)[2:end-1]
            maximum(abs.(re_diff)) < 1e-3 || return false
            maximum(abs.(im_diff)) < 1e-3 || return false
        end
    end

    return true
end

function sharpboundaryreflection_deriv()
    eas, freq, ground, bfield, electrons = derivative_scenario()

    M = LWMS.susceptibility(70e3, freq, bfield, electrons)

    Ros = SMatrix{2,2,ComplexF64,4}[]
    dRs = SMatrix{2,2,ComplexF64,4}[]
    for ea in eas
        Ro = LWMS.sharpboundaryreflection(ea, M)
        RdR = LWMS.sharpboundaryreflection(ea, M, LWMS.Derivative_dθ())

        Ro ≈ RdR[SVector(1,2),:] || return false

        push!(Ros, Ro)
        push!(dRs, RdR[SVector(3,4),:])
    end

    re_ros = Tuple.(real.(Ros))
    im_ros = Tuple.(imag.(Ros))

    re_drs = Tuple.(real.(dRs))
    im_drs = Tuple.(imag.(dRs))

    for j = 1:4
        re_diff = centereddiff(getfield.(re_ros, j), Δθ) - getfield.(re_drs, j)[2:end-1]
        im_diff = centereddiff(getfield.(im_ros, j), Δθ) - getfield.(im_drs, j)[2:end-1]
        maximum(abs.(re_diff)) < 1e-5 || return false
        maximum(abs.(im_diff)) < 1e-5 || return false
    end

    return true
end

# NOTE: This function takes a relatively long time because it has to do 2*length(eas) integrations
function integratedreflection_deriv()
    eas, freq, ground, bfield, electrons = derivative_scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    Ros = SMatrix{2,2,ComplexF64,4}[]
    dRs = SMatrix{2,2,ComplexF64,4}[]
    @info "    Integrating dR/dz twice for each eigenangle. This may take a minute."
    for ea in eas
        Ro = LWMS.integratedreflection(ea, freq, waveguide)
        RdR = LWMS.integratedreflection(ea, freq, waveguide, LWMS.Derivative_dθ())

        isapprox(Ro, RdR[SVector(1,2),:], rtol=1e-3) || return false

        push!(Ros, Ro)
        push!(dRs, RdR[SVector(3,4),:])
    end
    @info "    Integrations of dR/dz complete."

    re_ros = Tuple.(real.(Ros))
    im_ros = Tuple.(imag.(Ros))

    re_drs = Tuple.(real.(dRs))
    im_drs = Tuple.(imag.(dRs))

    for j = 1:4
        # XXX: Tolerance is so low because R we might be passing through a branch point and
        # R is unreasonably large (and probably) not defined
        # TODO: Handle branch point calculation of R
        isapprox(centereddiff(getfield.(re_ros, j), Δθ), getfield.(re_drs, j)[2:end-1], rtol=0.1) || return false
        isapprox(centereddiff(getfield.(im_ros, j), Δθ), getfield.(im_drs, j)[2:end-1], rtol=0.1) || return false
    end

    return true
end

function fresnelreflection_deriv()
    eas, freq, ground, bfield, electrons = derivative_scenario()

    Rgos = SMatrix{2,2,ComplexF64,4}[]
    dRgs = SMatrix{2,2,ComplexF64,4}[]
    for ea in eas
        Rgo = LWMS.fresnelreflection(ea, ground, freq)
        Rg, dRg = LWMS.fresnelreflection(ea, ground, freq, LWMS.Derivative_dθ())

        Rgo ≈ Rg || return false

        push!(Rgos, Rgo)
        push!(dRgs, dRg)
    end

    re_rgos = Tuple.(real.(Rgos))
    im_rgos = Tuple.(imag.(Rgos))

    re_drgs = Tuple.(real.(dRgs))
    im_drgs = Tuple.(imag.(dRgs))

    for j = 1:4
        re_diff = centereddiff(getfield.(re_rgos, j), Δθ) - getfield.(re_drgs, j)[2:end-1]
        im_diff = centereddiff(getfield.(im_rgos, j), Δθ) - getfield.(im_drgs, j)[2:end-1]
        maximum(abs.(re_diff)) < 1e-4 || return false
        maximum(abs.(im_diff)) < 1e-4 || return false
    end

    return true
end

function modalequation_deriv()
    eas, freq, ground, bfield, electrons = derivative_scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    Fs = ComplexF64[]
    dFs = ComplexF64[]
    for ea in eas
        RdR = LWMS.integratedreflection(ea, freq, waveguide, LWMS.Derivative_dθ())
        Rg, dRg = LWMS.fresnelreflection(ea, ground, freq, LWMS.Derivative_dθ())

        R = RdR[SVector(1,2),:]
        dR = RdR[SVector(3,4),:]

        F = LWMS.modalequation(R, Rg)
        dF = LWMS.modalequationdθ(R, dR, Rg, dRg)

        push!(Fs, F)
        push!(dFs, dF)
    end
    # XXX: Suffering from same branch point (?) issue as `eval_integratedreflection`
    isapprox(centereddiff(real.(Fs), Δθ), real.(dFs)[2:end-1], rtol=0.05) || return false
    isapprox(centereddiff(imag.(Fs), Δθ), imag.(dFs)[2:end-1], rtol=0.05) || return false

    return true
end

function resonantmodalequation_deriv()
    eas, freq, ground, bfield, electrons = derivative_scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    # known solution
    ea = EigenAngle(0.8654508431405412 - 0.026442282713736928im)
    dFdθ, R, Rg = LWMS.solvemodalequationdθ(ea, freq, waveguide)

    # If dθ is too small, the finitedF blows up
    dθ = deg2rad(1e-2)
    Fm = LWMS.solvemodalequation(EigenAngle(ea.θ-dθ/2), freq, waveguide)
    Fp = LWMS.solvemodalequation(EigenAngle(ea.θ+dθ/2), freq, waveguide)

    finitedF = (Fp - Fm)/dθ

    return isapprox(dFdθ, finitedF, rtol=1e-2)
end

########
# No B-field
########



########
# Vertical B-field
########
# @testset "Vertical B-field" begin
#

function scenario()
    ea = EigenAngle(deg2rad(complex(85.0, -1.0)))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    # `tx` isn't just bits
    ground = Ground(15, 0.001)
    bfield = BField(50e-6, deg2rad(90), 0)

    electrons = Species(qₑ, mₑ,
                        z -> waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
                        z -> electroncollisionfrequency(z, LWMS.CURVATURE_HEIGHT))

    return ea, tx, ground, bfield, electrons
end

function test_wmatrix()
    # Make sure "analytical" solution for W matches numerical
    ea, tx, ground, bfield, electrons = scenario()

    C = ea.cosθ
    L = @SMatrix [C 0 -C 0;
                  0 -1 0 -1;
                  0 -C 0 C;
                  1 0 1 0]

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)

    return [W[1] W[3]; W[2] W[4]] ≈ 2*(L\T)*L
end

function test_bookerquartic()
    ea, tx, ground, bfield, electrons = scenario()

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)
    LWMS.bookerquartic!(ea, M)

    q, B = LWMS.BOOKER_QUARTIC_ROOTS, LWMS.BOOKER_QUARTIC_COEFFS

    S = ea.sinθ
    S², C² = ea.sin²θ, ea.cos²θ

    # Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
    for i in eachindex(q)
        G = [1-q[i]^2+M[1,1] M[1,2] S*q[i]+M[1,3];
             M[2,1] 1-q[i]^2-S²+M[2,2] M[2,3];
             S*q[i]+M[3,1] M[3,2] C²+M[3,3]]
        isapprox(det(G), 0, atol=sqrt(eps())) || return false
    end

    # Confirm Booker quartic is directly satisfied
    for i in eachindex(q)
        booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
        isapprox(booker, 0, atol=sqrt(eps())) || return false
    end

    return true
end

function test_sharplybounded()
    ea, tx, ground, bfield, electrons = scenario()

    M = LWMS.susceptibility(95e3, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)
    LWMS.bookerquartic!(ea, M)

    initR = LWMS.sharpboundaryreflection(ea, M)

    initR[1,2] ≈ initR[2,1] || return false

    function iterativesharpR!(f, R, W)
        f .= W[2] + W[4]*R - R*W[1] - R*W[3]*R
    end
    initR0 = complex([1 0.1; 0.1 1])
    res = nlsolve((f,x)->iterativesharpR!(f,x,W), initR0)

    return initR ≈ res.zero
end

# BUG: This isn't working? Possibly foundational math error
function numericalsharpR()
    ea, tx, ground, bfield, electrons = scenario()

    M = LWMS.susceptibility(95e3, tx.frequency, bfield, electrons)
    LWMS.bookerquartic!(ea, M)
    q, B = LWMS.BOOKER_QUARTIC_ROOTS, LWMS.BOOKER_QUARTIC_COEFFS

    sort!(q, by=LWMS.upgoing)

    S, C = ea.sinθ, ea.cosθ
    S², C² = ea.sin²θ, ea.cos²θ
    esum = zeros(ComplexF64,4)
    for i in 1:2
        Γ = [0 -q[i] 0;
             q[i] 0 -S;
             0 S 0]
        G = Γ^2 + I + M
        E = nullspace(G)
        ℋ = Γ*E  # Γℋ == -(I + M)*E

        esum += [E[1]; E[2]; ℋ[1]; ℋ[2]]
    end
    A = [C 0 -C 0; 0 1 0 1; 0 -C 0 C; 1 0 1 0]
    e = A\esum

    R = [e[3]/e[1] e[4]/e[1];
         e[3]/e[2] e[4]/e[2]]

    initR = LWMS.sharpboundaryreflection(ea, M)

    initR ≈ R
end

function verticalreflection()
    ea, tx, ground, bfield, electrons = scenario()

    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    R = LWMS.integratedreflection(ea, tx.frequency, waveguide)

    return R[1,2] ≈ R[2,1]
end

function pecground()
    pec_ground = LWMS.Ground(1, 1e12)
    vertical_ea = LWMS.EigenAngle(0)  # not necessary for PEC test?

    return abs.(LWMS.fresnelreflection(vertical_ea, pec_ground, Frequency(24e3))) ≈ I
end

function modalequation()
    ea, tx, ground, bfield, electrons = scenario()

    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    # TODO: check this - it's not a great match to zero
    # known solution w/ tolerance=1e-9
    ea = EigenAngle(1.305895494889554 - 0.030739068016986393im)
    f = LWMS.solvemodalequation(ea, tx.frequency, waveguide)

    return isapprox(f, complex(0), atol=0.01)
end

function modefinder()
    ea, tx, ground, bfield, electrons = scenario()

    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    # TODO: Make a `triangulardomain` for this problem o avoid the low right
    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8

    modes = LWMS.findmodes(origcoords, tx.frequency, waveguide, tolerance)

    for m in modes
        f = LWMS.solvemodalequation(EigenAngle(m), tx.frequency, waveguide)
        isapprox(f, complex(0), atol=1e-4) || return false
    end
    return true
end

@testset "modefinder.jl" begin
    @info "Testing modefinder"

    @test test_wmatrix()
    @test test_bookerquartic()
    @test test_sharplybounded()
    @test_skip numericalsharpR()
    @test verticalreflection()
    @test pecground()
    @test modalequation()
    @test modefinder()

    # TODO: Off-diagonal terms should be 0 with no B field
    @test_skip +(M[1,2], M[1,3], M[2,1], M[2,3], M[3,1], M[3,2]) == 0

    @testset "Derivatives" begin
        @info "  Derivatives..."

        @test wmatrix_deriv()
        @test sharpboundaryreflection_deriv()
        @test fresnelreflection_deriv()
        @test integratedreflection_deriv()
        @test modalequation_deriv()
        @test resonantmodalequation_deriv()
    end
end
