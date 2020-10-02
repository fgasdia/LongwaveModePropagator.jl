err_func(a,b) = maximum(abs.(a-b))

function derivative_scenario()
    θs = [complex(r,i) for r = range(deg2rad(30), deg2rad(89), length=100) for i = range(deg2rad(-30), deg2rad(0), length=100)]

    freq = Frequency(24e3)
    bfield = BField(50e-6, π/2, 0)
    ground = Ground(15, 0.001)
    electrons = Species(QE, ME,
                        z -> waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
                        z -> electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT))

    return θs, freq, ground, bfield, electrons
end

function wmatrix_deriv()
    θs, freq, ground, bfield, electrons = derivative_scenario()

    M = LWMS.susceptibility(70e3, freq, bfield, electrons)

    for i = 1:4
        for j = 1:4
            Wfcn(θ) = (ea = EigenAngle(θ); T = LWMS.tmatrix(ea, M); LWMS.wmatrix(ea, T)[i][j])
            dWref = FiniteDiff.finite_difference_derivative(Wfcn, θs, Val{:central})
            dW(θ) = (ea = EigenAngle(θ); T = LWMS.tmatrix(ea, M); LWMS.dwmatrixdθ(ea, M, T)[i][j])

            err_func(dW.(θs), dWref) < 1e-3 || return false
        end
    end

    return true
end

function sharpboundaryreflection_deriv()
    θs, freq, ground, bfield, electrons = derivative_scenario()

    M = LWMS.susceptibility(70e3, freq, bfield, electrons)

    for i = 1:4
        Rref(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M)[i])
        dRref = FiniteDiff.finite_difference_derivative(Rref, θs, Val{:central})
        R(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M, LWMS.Derivative_dθ())[SVector(1,2),:][i])
        dR(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M, LWMS.Derivative_dθ())[SVector(3,4),:][i])

        Rref.(θs) ≈ R.(θs) || return false
        err_func(dR.(θs), dRref) < 1e-4 || return false
    end

    return true
end

function integratedreflection_deriv()
    θs, freq, ground, bfield, electrons = derivative_scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(alt) = LWMS.susceptibility(alt, freq, bfield, electrons)

    @info "    Integrating dR/dz for many eigenangles. This may take a minute."

    Rref(θ) = (ea = EigenAngle(θ); LWMS.integratedreflection(ea, freq, waveguide, Mfcn)[1])
    dRref = FiniteDiff.finite_difference_derivative(Rref, θs, Val{:central})
    R(θ) = (ea = EigenAngle(θ); LWMS.integratedreflection(ea, freq, waveguide, Mfcn, LWMS.Derivative_dθ())[SVector(1,2),:][1])
    dR(θ) = (ea = EigenAngle(θ); LWMS.integratedreflection(ea, freq, waveguide, Mfcn, LWMS.Derivative_dθ())[SVector(3,4),:][1])

    # Rref.(θs) !≈ R.(θs) because of difference in integration atol. For a small number of
    # `R` values that are quite large (thousands), the difference in tolerance can be significant
    Rdiff = Rref.(θs) .- R.(θs)
    quantile(real(Rdiff), 0.99) < 0.01 || return false
    quantile(imag(Rdiff), 0.99) < 0.01 || return false

    # Similar issue with dR...
    dRθs = dR.(θs)
    dRdiff = dRθs .- dRref
    quantile(abs.(dRdiff./dRref), 0.99) < 1 || return false

    return true
end

function fresnelreflection_deriv()
    θs, freq, ground, bfield, electrons = derivative_scenario()

    for i = 1:4
        Rgref(θ) = (ea = EigenAngle(θ); LWMS.fresnelreflection(ea, ground, freq)[1])
        dRgref = FiniteDiff.finite_difference_derivative(Rgref, θs, Val{:central})
        Rg(θ) = (ea = EigenAngle(θ); LWMS.fresnelreflection(ea, ground, freq, LWMS.Derivative_dθ())[1][1])
        dRg(θ) = (ea = EigenAngle(θ); LWMS.fresnelreflection(ea, ground, freq, LWMS.Derivative_dθ())[2][1])

        Rgref.(θs) ≈ Rg.(θs) || return false
        err_func(dRg.(θs), dRgref) < 1e-6 || return false
    end

    return true
end

function modalequation_deriv()
    θs, freq, ground, bfield, electrons = derivative_scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(alt) = LWMS.susceptibility(alt, freq, bfield, electrons)

    function Fref(θ)
        ea = EigenAngle(θ)
        R = LWMS.integratedreflection(ea, freq, waveguide, Mfcn)
        Rg = LWMS.fresnelreflection(ea, ground, freq)

        return LWMS.modalequation(R, Rg)
    end

    function dF(θ)
        ea = EigenAngle(θ)
        RdR = LWMS.integratedreflection(ea, freq, waveguide, Mfcn, LWMS.Derivative_dθ())
        Rg, dRg = LWMS.fresnelreflection(ea, ground, freq, LWMS.Derivative_dθ())

        R = RdR[SVector(1,2),:]
        dR = RdR[SVector(3,4),:]

        return LWMS.modalequationdθ(R, dR, Rg, dRg)
    end

    @info "    Integrating dR/dz for many eigenangles. This may take a minute."

    Fθ = Fref.(θs)
    dFref = FiniteDiff.finite_difference_derivative(Fref, θs, Val{:forward}, eltype(θs), Fθ)

    dFθ = dF.(θs)
    dFdiff = dFθ .- dFref
    quantile(abs.(dFdiff./dFref), 0.99) < 2 || return false

    return true
end

function resonantmodalequation_deriv()
    θs, freq, ground, bfield, electrons = derivative_scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    Mfcn(alt) = LWMS.susceptibility(alt, freq, bfield, electrons)
    modeequation = LWMS.PhysicalModeEquation(freq, waveguide, Mfcn)

    # known solution
    θ = 1.4152764714690873 - 0.01755942938376613im
    ea = EigenAngle(θ)
    dFdθ, R, Rg = LWMS.solvemodalequation(ea, modeequation, LWMS.Derivative_dθ())

    Fref(θ) = (ea = EigenAngle(θ); LWMS.solvemodalequation(ea, modeequation))
    dFref = FiniteDiff.finite_difference_derivative(Fref, θ, Val{:central})

    return isapprox(dFdθ, dFref, rtol=1e-1)
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
    bfield = BField(50e-6, π/2, 0)

    electrons = Species(QE, ME,
                        z -> waitprofile(z, 75, 0.32, cutoff_low=LWMS.CURVATURE_HEIGHT),
                        z -> electroncollisionfrequency(z, cutoff_low=LWMS.CURVATURE_HEIGHT))

    return ea, tx, ground, bfield, electrons
end

function test_wmatrix()
    # Check that "analytical" solution for W matches numerical
    ea, tx, ground, bfield, electrons = scenario()

    C = ea.cosθ
    L = [C 0 -C 0;
         0 -1 0 -1;
         0 -C 0 C;
         1 0 1 0]

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)

    return [W[1] W[3]; W[2] W[4]] ≈ 2*(L\T)*L
end

function test_sharplybounded()
    ea, tx, ground, bfield, electrons = scenario()

    M = LWMS.susceptibility(95e3, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)

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
    q, B = LWMS.bookerquartic(ea, M)

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
        H = Γ*E  # Γℋ == -(I + M)*E

        esum += [E[1]; E[2]; H[1]; H[2]]
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

    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
    R = LWMS.integratedreflection(ea, tx.frequency, waveguide, Mfcn)

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
    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

    # known solution
    ea = EigenAngle(1.4152764714690873 - 0.01755942938376613im)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
    f = LWMS.solvemodalequation(ea, modeequation)

    return LWMS.isroot(f)
end

function modefinder()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, electrons)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)

    # Δr from 0.5->0.25 => time from 3.8->5.3 sec
    # tolerance from 1e-8->1e-7 => time from 5.3->4.6 sec
    origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = LWMS.GRPFParams(est_num_nodes, 1e-6, true)

    modes = LWMS.findmodes(origcoords, grpfparams, modeequation)
    
    for m in modes
        f = LWMS.solvemodalequation(m, modeequation)
        # println(f)
        LWMS.isroot(f) || return false
    end
    return true
end













# TEMP function for testing
function manyint(eas)
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
    for i in eachindex(eas)
        LWMS.integratedreflection(ea, tx.frequency, waveguide, Mfcn)
    end
end

function manysolves(eas)
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
    # Mfcn = LWMS.susceptibilityinterpolator(tx.frequency, bfield, electrons)

    # eas = EigenAngle.(rand(ComplexF64, 1000))

    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
    for i in eachindex(eas)
        LWMS.solvemodalequation(eas[i], modeequation)
    end
end

function findmodes()
    ea, tx, ground, bfield, electrons = scenario()
    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

    # Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, electrons)
    Mfcn = LWMS.susceptibilityinterpolator(tx.frequency, bfield, electrons)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)

    # Δr from 0.5->0.25 => time from 3.8->5.3 sec
    # tolerance from 1e-8->1e-7 => time from 5.3->4.6 sec
    origcoords = LWMS.defaultcoordinates(tx.frequency.f)
    est_num_nodes = ceil(Int, length(origcoords)*1.5)
    grpfparams = LWMS.GRPFParams(est_num_nodes, 1e-5, true)

    f(θ) = LWMS.solvemodalequation(EigenAngle(θ), modeequation)
    roots, poles, quads, diffs, tess, g2f = LWMS.grpf(f, origcoords, LWMS.PlotData(), grpfparams)

    # Ensure roots are valid solutions to the modal equation
    LWMS.filterroots!(roots, tx.frequency, waveguide)

    # # Remove any redundant modes
    # if tolerance is 1e-8, this rounds to 7 decimal places
    sort!(roots, by=reim)
    ndigits = round(Int, abs(log10(grpfparams.tolerance)+1), RoundDown)
    unique!(z->round(z, digits=ndigits), roots)

    # return EigenAngle.(roots)
    return roots, est_num_nodes, length(quads)
end


function customsolves(eas, tx, waveguide)
    @unpack bfield, species = waveguide
    frequency = tx.frequency
    Mfcn(z) = LWMS.susceptibility(z, tx.frequency, bfield, species)

    for ea in eas
        R = LWMS.integratedreflection(ea, frequency, waveguide, Mfcn)
        Rg = LWMS.fresnelreflection(ea, waveguide.ground, frequency)

        f = LWMS.modalequation(R, Rg)
    end
end

function test_triangulardomain()
    za = complex(60, -0.1)
    zb = complex(89, -0.1)
    zc = complex(60, -10.0)
    r = 1
    v = LWMS.triangulardomain(za, zb, zc, r)
    # plot(real(v),imag(v),"o")

    # Scale points for tesselation
    rmin, rmax = minimum(real, v), maximum(real, v)
    imin, imax = minimum(imag, v), maximum(imag, v)

    width = RootsAndPoles.MAXCOORD - RootsAndPoles.MINCOORD
    ra = width/(rmax-rmin)
    rb = RootsAndPoles.MAXCOORD - ra*rmax

    ia = width/(imax-imin)
    ib = RootsAndPoles.MAXCOORD - ia*imax

    nodes = [Point2D(real(coord), imag(coord)) for coord in v]
    tess = DelaunayTessellation(length(nodes))
    push!(tess, nodes)
end

@testset "modefinder.jl" begin
    @info "Testing modefinder"

    @test test_wmatrix()
    @test test_sharplybounded()
    @test_skip numericalsharpR()
    @test verticalreflection()
    @test pecground()
    @test modalequation()
    @test modefinder()

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
