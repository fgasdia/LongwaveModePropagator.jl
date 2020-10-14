function wmatrix_deriv(scenario)
    @unpack tx, bfield, species = scenario

    M = LWMS.susceptibility(70e3, tx.frequency, bfield, species)

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

function sharpboundaryreflection_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)

    for i = 1:4
        Rref(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M)[i])
        dRref = FiniteDiff.finite_difference_derivative(Rref, θs, Val{:central})
        R(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M, LWMS.Dθ())[SVector(1,2),:][i])
        dR(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M, LWMS.Dθ())[SVector(3,4),:][i])

        Rref.(θs) ≈ R.(θs) || return false
        err_func(dR.(θs), dRref) < 1e-6 || return false
    end

    return true
end

function sharpboundaryreflection_sheddy(ea, M, ::LWMS.Dθ)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ
    C2 = 2*C

    q, B = LWMS.bookerquartic(ea, M)

    # We choose the 2 roots corresponding to upward travelling waves as being
    # those that lie close to the positive real axis and negative imaginary axis
    LWMS.sortquarticroots!(q)

    #==
    Γ = [0 -q 0;
         q 0 -S;
         0 S 0]
    G = Γ² + I + M
    G = [1-q[i]^2+M[1,1] M[1,2] S*q[i]+M[1,3];
         M[2,1] 1-q[i]^2-S²+M[2,2] M[2,3];
         S*q[i]+M[3,1] M[3,2] C²+M[3,3]]
    ==#

    # NOTE: Sheddy specifies roots 2, 1 rather than the order used by Pitteway

    # Precompute
    q₁S = q[1]*S
    q₂S = q[2]*S

    M11p1 = 1 + M[1,1]

    # Constant entries of dispersion matrix `G`
    G12 = M[1,2]
    G32 = M[3,2]
    G33 = C² + M[3,3]

    G12G33 = G12*G33

    # Values for solutions of the Booker quartic corresponding to upgoing waves
    G11₁ = M11p1 - q[1]^2
    G13₁ = M[1,3] + q₁S
    G31₁ = M[3,1] + q₁S

    Δ₁ = G11₁*G33 - G13₁*G31₁
    Δ₁⁻¹ = 1/Δ₁
    P₁ = (-G12G33 + G13₁*G32)*Δ₁⁻¹
    T₁ = q[1]*P₁ - S*(-G11₁*G32 + G12*G31₁)*Δ₁⁻¹
    T₁C = T₁*C

    G11₂ = M11p1 - q[2]^2
    G13₂ = M[1,3] + q₂S
    G31₂ = M[3,1] + q₂S

    Δ₂ = G11₂*G33 - G13₂*G31₂
    Δ₂⁻¹ = 1/Δ₂
    P₂ = (-G12G33 + G13₂*G32)*Δ₂⁻¹
    T₂ = q[2]*P₂ - S*(-G11₂*G32 + G12*G31₂)*Δ₂⁻¹
    T₂C = T₂*C

    Δ = (T₁C + P₁)*(C + q[2]) - (T₂C + P₂)*(C + q[1])
    Δ⁻¹ = 1/Δ

    M13pM31 = M[1,3] + M[3,1]

    # Additional calculations required for dR/dθ
    dS = C
    dC = -S
    dC² = -S*C2
    dB3 = dS*M13pM31
    dB2 = -dC²*(2 + M[1,1] + M[3,3])
    dB1 = dS/S*B[2] - S*dC²*M13pM31
    dB0 = dC²*(2*C²*(1 + M[1,1]) + M[3,3] + M[2,2] + M[1,1]*(M[3,3] + M[2,2]) -
            M[1,3]*M[3,1] - M[1,2]*M[2,1])

    dq_1 = -(((dB3*q[1] + dB2)*q[1] + dB1)*q[1] + dB0) /
            (((4*B[5]*q[1] + 3*B[4])*q[1] + 2*B[3])*q[1] + B[2])
    dq_2 = -(((dB3*q[2] + dB2)*q[2] + dB1)*q[2] + dB0) /
            (((4*B[5]*q[2] + 3*B[4])*q[2] + 2*B[3])*q[2] + B[2])

    dG33 = dC²

    dG11₁ = -2*q[1]*dq_1
    dG13₁ = dq_1*S + q[1]*dS
    dG31₁ = dG13₁  # dq_1*S + q[1]*dS

    dΔ₁ = dG11₁*G33 + G11₁*dG33 - dG13₁*G31₁ - G13₁*dG31₁
    dΔ₁⁻¹ = -dΔ₁/Δ₁^2

    dP₁ = (-G12*dG33 + dG13₁*G32)*Δ₁⁻¹ + (G13₁*G32 - G12*G33)*dΔ₁⁻¹
    dT₁ = dq_1*P₁ + q[1]*dP₁ -
        dS*(-G11₁*G32 + G12*G31₁)*Δ₁⁻¹ -
        S*(-dG11₁*G32 + G12*dG31₁)*Δ₁⁻¹ -
        S*(-G11₁*G32 + G12*G31₁)*dΔ₁⁻¹

    dG11₂ = -2*q[2]*dq_2
    dG13₂ = dq_2*S + q[2]*dS
    dG31₂ = dG13₂  # dq_2*S + q[2]*dS

    dΔ₂ = dG11₂*G33 + G11₂*dG33 - dG13₂*G31₂ - G13₂*dG31₂
    dΔ₂⁻¹ = -dΔ₂/Δ₂^2

    dP₂ = (-G12*dG33 + dG13₂*G32)*Δ₂⁻¹ + (G13₂*G32 - G12*G33)*dΔ₂⁻¹
    dT₂ = dq_2*P₂ + q[2]*dP₂ -
        dS*(-G11₂*G32 + G12*G31₂)*Δ₂⁻¹ -
        S*(-dG11₂*G32 + G12*dG31₂)*Δ₂⁻¹ -
        S*(-G11₂*G32 + G12*G31₂)*dΔ₂⁻¹

    dΔ = dT₁*C² + T₁*dC² + dT₁*C*q[2] + T₁*dC*q[2] + T₁C*dq_2 + dP₁*C + P₁*dC +
            dP₁*q[2] + P₁*dq_2 - (dT₂*C² + T₂*dC²) -
            (dT₂*C*q[1] + T₂*dC*q[1] + T₂C*dq_1) -
            (dP₂*C + P₂*dC) - (dP₂*q[1] + P₂*dq_1)
    dΔ⁻¹ = -dΔ/Δ^2

    # R
    R11 = ((T₁C - P₁)*(C + q[2]) - (T₂C - P₂)*(C + q[1]))*Δ⁻¹  # ∥R∥
    R22 = ((T₁C + P₁)*(C - q[2]) - (T₂C + P₂)*(C - q[1]))*Δ⁻¹  # ⟂R⟂
    R12 = -C2*(T₁*P₂ - T₂*P₁)*Δ⁻¹  # ⟂R∥
    R21 = -C2*(q[1] - q[2])*Δ⁻¹  # ∥R⟂

    # dR
    dR11 = dΔ⁻¹*R11*Δ +
        Δ⁻¹*(C²*(dT₁ - dT₂) + dC²*(T₁ - T₂) + dC*(T₁*q[2] - P₁ - T₂*q[1] + P₂) +
            C*(dT₁*q[2] + T₁*dq_2 + dP₂ - dP₁ - dT₂*q[1] - T₂*dq_1) +
            dP₂*q[1] + P₂*dq_1 - dP₁*q[2] - P₁*dq_2)
    dR12 = dΔ⁻¹*R12*Δ + dC/C*R12 - C2*(dT₁*P₂ + T₁*dP₂ - dT₂*P₁ - T₂*dP₁)*Δ⁻¹
    dR21 = dΔ⁻¹*R21*Δ + dC/C*R21 - C2*(dq_1 - dq_2)*Δ⁻¹
    dR22 = dΔ⁻¹*R22*Δ +
            Δ⁻¹*(C²*(dT₁ - dT₂) + dC²*(T₁ - T₂) + dC*(T₂*q[1] - T₁*q[2] + P₁ - P₂) +
            C*(dP₁ - dT₁*q[2] + dT₂*q[1] - T₁*dq_2 + T₂*dq_1 - dP₂) +
            dP₂*q[1] + P₂*dq_1 - dP₁*q[2] - P₁*dq_2)

    #==
    [R11 R12;
     R21 R22;
     dR11 dR12;
     dR21 dR22]
    ==#
    return SMatrix{4,2}(R11, R21, dR11, dR21, R12, R22, dR12, dR22)
end

function sharpboundary_sheddy_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)

    sharpboundaryreflection_sheddy(ea, M, LWMS.Dθ()) ≈ LWMS.sharpboundaryreflection(ea, M, LWMS.Dθ())
end

function integratedreflection_deriv(scenario)
    @unpack tx, ground, bfield, species = scenario
    freq = tx.frequency

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(freq, waveguide)

    @info "    Integrating dR/dz for many eigenangles. This may take a minute."

    params = LWMS.IntegrationParams(1e-8, LWMS.BS5())
    Rref(θ,m) = (m = LWMS.setea(EigenAngle(θ), m); LWMS.integratedreflection(m, params=params))
    RdR(θ,m) = (m = LWMS.setea(EigenAngle(θ), m); LWMS.integratedreflection(m, LWMS.Dθ(), params=params))

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
        R = [v[i] for v in Rs]
        dR = [v[i] for v in dRs]
        Rr = [v[i] for v in Rrefs]
        dRr = [v[i] for v in dRrefs]

        # Rref.(θs) !≈ R.(θs) because of difference in integration tolerance.
        # For very large Rs, the difference can be significant, therefore we use rtol
        isapprox(R, Rr, rtol=1e-4) || return false
        isapprox(dR, dRr, rtol=1e-3) || return false
    end

    return true
end

function fresnelreflection_deriv(scenario)
    @unpack tx, ground = scenario
    freq = tx.frequency

    for i = 1:4
        Rgref(θ) = (ea = EigenAngle(θ); LWMS.fresnelreflection(ea, ground, freq)[i])
        dRgref = FiniteDiff.finite_difference_derivative(Rgref, θs, Val{:central})
        Rg(θ) = (ea = EigenAngle(θ); LWMS.fresnelreflection(ea, ground, freq, LWMS.Dθ())[1][i])
        dRg(θ) = (ea = EigenAngle(θ); LWMS.fresnelreflection(ea, ground, freq, LWMS.Dθ())[2][i])

        Rgref.(θs) ≈ Rg.(θs) || return false
        err_func(dRg.(θs), dRgref) < 1e-6 || return false
    end

    return true
end

function modalequation_deriv(scenario)
    @unpack ea, tx, ground, bfield, species = scenario
    freq = tx.frequency

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(ea, freq, waveguide)

    @info "    Integrating dR/dz for many eigenangles. This may take a minute."

    dFref = FiniteDiff.finite_difference_derivative(z->LWMS.solvemodalequation(z, modeequation),
        θs, Val{:central})

    dFs = Vector{ComplexF64}(undef, length(θs))
    Threads.@threads for i in eachindex(θs)
        dFdθ, R, Rg = LWMS.solvemodalequation(θs[i], modeequation, LWMS.Dθ())
        dFs[i] = dFdθ
    end

    return isapprox(dFs, dFref, rtol=1e-3)
end

########
# Non-derivative
########

function test_wmatrix(scenario)
    # Check that "analytical" solution for W matches numerical
    @unpack ea, tx, bfield, species = scenario

    C = ea.cosθ
    L = [C 0 -C 0;
         0 -1 0 -1;
         0 -C 0 C;
         1 0 1 0]

    M = LWMS.susceptibility(80e3, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)

    return [W[1] W[3]; W[2] W[4]] ≈ 2*(L\T)*L
end

function test_sharplybounded_vertical(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    initR = LWMS.sharpboundaryreflection(ea, M)

    return initR[1,2] ≈ initR[2,1]
end

function test_sharplybounded(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)

    initR = LWMS.sharpboundaryreflection(ea, M)

    function iterativesharpR!(f, R, W)
        f .= W[2] + W[4]*R - R*W[1] - R*W[3]*R
    end
    initR0 = complex([1 0.1; 0.1 -1])
    res = nlsolve((f,R)->iterativesharpR!(f,R,W), initR0)

    return initR ≈ res.zero
end

function verticalreflection(scenario)
    @unpack ea, tx, ground, bfield, species = scenario

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide)

    R = LWMS.integratedreflection(modeequation)

    return R[1,2] ≈ R[2,1]
end

function pecground()
    pec_ground = LWMS.Ground(1, 1e12)
    vertical_ea = LWMS.EigenAngle(0)  # not necessary for PEC test?

    return abs.(LWMS.fresnelreflection(vertical_ea, pec_ground, Frequency(24e3))) ≈ I
end

function modalequation(scenario)
    @unpack ea, tx, ground, bfield, species = scenario

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(ea, tx.frequency, waveguide)

    f = LWMS.solvemodalequation(modeequation)

    return LWMS.isroot(f)
end

function parallel_integration(scenario, eas)
    @unpack ea, tx, ground, bfield, species = scenario

    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide)

    # n = 1000
    # eas = EigenAngle.(complex.(π/2 .- rand(n), -rand(n)/10));
    Threads.@threads for i in eachindex(eas)
        modeequation = LWMS.setea(eas[i], modeequation)
        R = LWMS.integratedreflection(modeequation)
    end
end

function modefinder(scenario)
    @unpack tx, bfield, species, ground = scenario
    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide)

    origcoords = LWMS.defaultcoordinates(tx.frequency)
    # est_num_nodes = ceil(Int, length(origcoords)*1.5)
    # grpfparams = LWMS.GRPFParams(est_num_nodes, 1e-8, true)

    modes = LWMS.findmodes(modeequation, origcoords)

    for m in modes
        f = LWMS.solvemodalequation(m, modeequation)
        LWMS.isroot(f) || return false
    end
    return true
end

# Is it worth refining a low tolerance GRPF solution with Roots? probably not
# but best method is `Roots.muller`
# using Roots
#
# function refineroots(root, scenario)
#     #verticalB_scenario
#     @unpack tx, bfield, species, ground = scenario
#     waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
#
#     Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, species)
#     modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
#     f(θ) = LWMS.solvemodalequation(EigenAngle(θ), modeequation)
#     # df(θ) = LWMS.solvemodalequation(EigenAngle(θ), modeequation, LWMS.Dθ())[1]
#
#     # D(f) = x->ForwardDiff.derivative(f, float(x))
#     # Roots.newton(f, D(f), root, atol=1e-6)
#     # find_zero(f, root)
#     Roots.muller(f, root, xrtol=1e-8)
#     # Roots.secant_method(f, root, atol=1e-6)
#     # Roots.find_zero(f, root, Roots.Order1(), atol=1e-6, verbose=true)
# end
#
# function evalroot(root, scenario)
#     @unpack tx, bfield, species, ground = scenario
#     waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
#
#     Mfcn(alt) = LWMS.susceptibility(alt, tx.frequency, bfield, species)
#     modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide, Mfcn)
#     LWMS.solvemodalequation(EigenAngle(root), modeequation)
# end

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

    @test pecground()

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        @test test_wmatrix(scn)
        @test test_sharplybounded(scn)
        @test_broken numericalsharpR(scn)
        @test modefinder(scn)
    end
    for scn in (resonant_scenario, )
        @test modalequation(scn)
    end
    for scn in (verticalB_scenario, )
        @test test_sharplybounded_vertical(scn)
        @test verticalreflection(scn)
    end

    @testset "Derivatives" begin
        @info "  Derivatives..."
        for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
            @test wmatrix_deriv(scn)
            @test sharpboundaryreflection_deriv(scn)
            @test sharpboundary_sheddy_deriv(scn)
            @test fresnelreflection_deriv(scn)
            @test integratedreflection_deriv(scn)
            @test modalequation_deriv(scn)
        end
    end
end
