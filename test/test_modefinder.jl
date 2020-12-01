function wmatrix_deriv(scenario)
    @unpack tx, bfield, species = scenario

    M = LWMS.susceptibility(70e3, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Wfcn(θ) = (ea = EigenAngle(θ); T = LWMS.tmatrix(ea, M); LWMS.wmatrix(ea, T)[i][j])
            dWref = FiniteDiff.finite_difference_derivative(Wfcn, θs, Val{:central})
            dW(θ) = (ea = EigenAngle(θ); T = LWMS.tmatrix(ea, M); LWMS.wmatrix(ea, M, T, LWMS.Dθ())[i][j])

            err_func(dW.(θs), dWref) < 1e-5 || return false
        end
    end

    return true
end

function sharpboundaryreflection_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(110e3, tx.frequency, bfield, species)

    for i = (1,4)  # TODO: off-diagonal terms are very sensitive, hard to test
        Rref(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M)[i])
        dRref = FiniteDiff.finite_difference_derivative(Rref, θs, Val{:central})
        R(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M, LWMS.Dθ())[SVector(1,2),:][i])
        dR(θ) = (ea = EigenAngle(θ); LWMS.sharpboundaryreflection(ea, M, LWMS.Dθ())[SVector(3,4),:][i])

        Rref.(θs) ≈ R.(θs) || return false

        # last θs has a mismatch
        isapprox(dRref[1:end-1], dR.(θs)[1:end-1], rtol=1e-6) || return false
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

    M = LWMS.susceptibility(110e3, tx.frequency, bfield, species)

    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)

    ref = sharpboundaryreflection_sheddy(ea, M, LWMS.Dθ())
    test = LWMS.sharpboundaryreflection(ea, M, LWMS.Dθ())

    # Off-diagonal terms of W are very small <1e-20 and off-diagonals of R are very sensitive
    diag(ref[1:2,:]) ≈ diag(test[1:2,:]) || return false
    isapprox(diag(ref[3:4,:]), diag(test[3:4,:]), rtol=1e-6) || return false

    return true
end

function sharpboundaryreflection_wavefields(ea::EigenAngle, M, ::LWMS.Dθ)
    S, C = ea.sinθ, ea.cosθ

    T = LWMS.tmatrix(ea, M)
    dT = LWMS.tmatrix(ea, M, LWMS.Dθ())
    q, B = LWMS.bookerquartic(T)
    LWMS.sortquarticroots!(q)
    dq = LWMS.bookerquartic(ea, M, q, B, LWMS.Dθ())

    # Compute fields using the ratios from [^Budden1988] pg 190
    a1 = -(T[1,1] + T[4,4])
    a2 = T[1,1]*T[4,4] - T[1,4]*T[4,1]
    a3 = T[1,2]
    a4 = T[1,4]*T[4,2] - T[1,2]*T[4,4]
    a5 = T[4,2]
    a6 = T[1,2]*T[4,1] - T[1,1]*T[4,2]

    # Some `dT` components are 0
    da1 = -(dT[1,1] + dT[4,4])
    da2 = T[1,1]*dT[4,4] + dT[1,1]*T[4,4] - dT[1,4]*T[4,1]
    da3 = dT[1,2]
    da4 = dT[1,4]*T[4,2] - T[1,2]*dT[4,4] - dT[1,2]*T[4,4]
    # da5 = 0
    da6 = dT[1,2]*T[4,1] - dT[1,1]*T[4,2]

    e = MMatrix{4,2,eltype(T),8}(undef)
    de = MMatrix{4,2,eltype(T),8}(undef)
    @inbounds for i = 1:2
        A = q[i]^2 + a1*q[i] + a2
        dA = 2*q[i]*dq[i] + a1*dq[i] + q[i]*da1 + da2

        e[1,i] = a3*q[i] + a4
        e[2,i] = A
        e[3,i] = q[i]*A
        e[4,i] = a5*q[i] + a6

        de[1,i] = a3*dq[i] + da3*q[i] + da4
        de[2,i] = dA
        de[3,i] = dq[i]*A + q[i]*dA
        de[4,i] = a5*dq[i] + da6  # + da5*q[i] = 0
    end

    d = SMatrix{2,2}(C*e[4,1]-e[1,1], -C*e[2,1]+e[3,1], C*e[4,2]-e[1,2], -C*e[2,2]+e[3,2])
    dd = SMatrix{2,2}(-S*e[4,1] + C*de[4,1] - de[1,1], S*e[2,1] - C*de[2,1] + de[3,1],
                      -S*e[4,2] + C*de[4,2] - de[1,2], S*e[2,2] - C*de[2,2] + de[3,2])
    u = SMatrix{2,2}(C*e[4,1]+e[1,1], -C*e[2,1]-e[3,1], C*e[4,2]+e[1,2], -C*e[2,2]-e[3,2])
    du = SMatrix{2,2}(-S*e[4,1] + C*de[4,1] + de[1,1], S*e[2,1] - C*de[2,1] - de[3,1],
                      -S*e[4,2] + C*de[4,2] + de[1,2], S*e[2,2] - C*de[2,2] - de[3,2])

    R = d/u
    dR = dd/u + d*(-u\du/u)  # BUG? dR[2,2] seems to be wrong...

    #==
    [R11 R12;
     R21 R22;
     dR11 dR12;
     dR21 dR22]
    ==#
    return vcat(R, dR)
end

function sharpboundary_wavefields_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LWMS.susceptibility(110e3, tx.frequency, bfield, species)

    T = LWMS.tmatrix(ea, M)
    W = LWMS.wmatrix(ea, T)

    ref = sharpboundaryreflection_wavefields(ea, M, LWMS.Dθ())
    test = LWMS.sharpboundaryreflection(ea, M, LWMS.Dθ())

    # Off-diagonal terms of W are very small <1e-20 and off-diagonals of R are very sensitive
    diag(ref[1:2,:]) ≈ diag(test[1:2,:]) || return false

    return true
end

function integratedreflection_deriv(scenario)
    @unpack tx, ground, bfield, species = scenario
    freq = tx.frequency

    # params = LWMSParams(integrationparams=IntegrationParams(1e-7, LWMS.RK4(), false))
    params = LWMSParams()
    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)
    modeequation = LWMS.PhysicalModeEquation(freq, waveguide)

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
        R = [v[i] for v in Rs[1:end-1]]
        dR = [v[i] for v in dRs[1:end-1]]
        Rr = [v[i] for v in Rrefs[1:end-1]]
        dRr = [v[i] for v in dRrefs[1:end-1]]

        isapprox(R, Rr, rtol=1e-4) || return false
        isapprox(dR, dRr, rtol=1e-3, atol=1e-3) || return false
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

    dFref = FiniteDiff.finite_difference_derivative(z->LWMS.solvemodalequation(z, modeequation),
        θs, Val{:central})

    dFs = Vector{ComplexF64}(undef, length(θs))
    Threads.@threads for i in eachindex(θs)
        dFdθ, R, Rg = LWMS.solvemodalequation(θs[i], modeequation, LWMS.Dθ())
        dFs[i] = dFdθ
    end

    return isapprox(dFs, dFref, rtol=1e-4)
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

    M = LWMS.susceptibility(110e3, tx.frequency, bfield, species)
    initR = LWMS.sharpboundaryreflection(ea, M)

    return initR[1,2] ≈ initR[2,1]
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
    vertical_ea = LWMS.EigenAngle(π/2)

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

    params = LWMSParams(grpfparams=LWMS.GRPFParams(100000, 1e-6, true))
    modes = LWMS.findmodes(modeequation, origcoords, params=params)

    for m in modes
        f = LWMS.solvemodalequation(m, modeequation, params=params)
        LWMS.isroot(f) || return false
    end
    return modes
    return true
end

function evalroot(root, scenario)
    @unpack tx, bfield, species, ground = scenario
    waveguide = LWMS.HomogeneousWaveguide(bfield, species, ground)

    modeequation = LWMS.PhysicalModeEquation(tx.frequency, waveguide)
    LWMS.solvemodalequation(EigenAngle(root), modeequation)
end


@testset "modefinder.jl" begin
    @info "Testing modefinder"

    @test pecground()

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        @test test_wmatrix(scn)
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
            @test sharpboundary_sheddy_deriv(scn)  # sensitive off-diag
            @test sharpboundary_wavefields_deriv(scn)
            @test fresnelreflection_deriv(scn)
            @test integratedreflection_deriv(scn) # again, a problem with off-diag
            @test modalequation_deriv(scn)
        end
    end
end
