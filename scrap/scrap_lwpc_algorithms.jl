
function what(v)
    angq = angle(v)
    imag(v) < 0 && (angq += 2π)
    angq < -π/4 && (angq += 2π)
    return angq
end

function what2(v)
    angq = rad2deg(angle(v))
    imag(v) < 0 && (angq += 360)
    angq < 315 && (angq += 360)
    return angq
end

function lwsort(v)
    angq = rad2deg(angle(v))
    angq < 0 && (angq += 360)
    angq < 135 && (angq += 360)
    return abs(angq-315)
end


function lwsortfull(v)
    q = copy(v)
    diffq = similar(q, Float64)

    l = 0
    for m = 2:4
        for n = m:4
            if imag(q[n]) <= 0
                l += 1
                q[n], q[m-1] = q[m-1], q[n]
            end
        end
    end
    if l != 2
        for n = 1:4
            diffq[n] = lwsort(q[n])
        end
        for nm = 2:4
            for n = nm:4
                if diffq[n] <= diffq[nm-1]
                    diffq[n], diffq[nm-1] = diffq[nm-1], diffq[n]
                    q[n], q[nm-1] = q[nm-1], q[n]
                end
            end
        end
    end
    return q
end

function wfsortfull(v)
    q = copy(v)
    diffq = similar(q, Float64)

    imag(q[1]) <= 0 ? l = 1 : l = 0

    m = l + 1
    while m < 4
        n = m + 1
        while n <= 4
            if imag(q[n]) <= 0
                l += 1
                q[n], q[m] = q[m], q[n]
                m += 1
            end
            n+= 1
        end
        m += 1
    end

    if l != 2
        m = 1
        while m <= 4
            diffq[m] = what2(q[m])
            m += 1
        end
        m = 2
        while m <= 4
            n = m
            while n <= 4
                if diffq[n] <= diffq[m-1]
                    diffq[n], diffq[m-1] = diffq[m-1], diffq[n]
                    q[n], q[m-1] = q[m-1], q[n]
                end
                n += 1
            end
            m += 1
        end
    end
    return q
end


function lwpcreflectioncoeffs(ea::EigenAngle, e1, e2)
    # From wf_r_mtrx.for

    C = ea.cosθ

    g12 = e1[1]*e2[2] - e2[1]*e1[2]
    g13 = e1[1]*e2[3] - e2[1]*e1[3]
    g14 = e1[1]*e2[4] - e2[1]*e1[4]
    g23 = e1[2]*e2[3] - e2[2]*e1[3]
    g24 = e1[2]*e2[4] - e2[2]*e1[4]
    g34 = e1[3]*e2[4] - e2[3]*e1[4]

    den = -g13 + C*(g34 - g12 + C*g24)

    d11 = g13 + C*(g34 + g12 + C*g24)
    d22 = g13 + C*(-g34 - g12 + C*g24)
    d12 = 2*C*g14
    d21 = 2*C*g23

    return SMatrix{2,2,eltype(den),4}(d11/den, d21/den, d12/den, d22/den)
end

vacuumR = vacuumreflectioncoeffs(ea, e1[end], e2[end])
@test vacuumR ≈ lwpcreflectioncoeffs(ea, e1[end], e2[end])

function lwpcscale(p1, p2)
    e1, e2 = MVector(p1), MVector(p2)

    # aterm → dot(e1, e1)
    aterm = 0
    for i = 1:4
        aterm += abs2(e1[i])
    end

    # term → dot(e1, e2)
    term = 0
    for i = 1:4
        term += conj(e1[i])*e2[i]
    end

    # term → dot(e1, e2)/dot(e1, e1)
    term /= aterm

    # e2 → e2 - dot(e1, e2)/dot(e1, e1)
    for i = 1:4
        e2[i] -= term*e1[i]
    end

    # bterm → dot(e2, e2)
    bterm = 0
    for i = 1:4
        bterm += abs2(e2[i])
    end

    # Normalize both vectors
    aterm = 1/sqrt(aterm)
    bterm = 1/sqrt(bterm)
    for i = 1:4
        e1[i] *= aterm
        e2[i] *= bterm
    end

    return SVector(e1), SVector(e2)
end



function sharpboundaryreflection_sheddy(ea, M, ::LMP.Dθ)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ
    C2 = 2*C

    q, B = LMP.bookerquartic(ea, M)

    # We choose the 2 roots corresponding to upward travelling waves as being
    # those that lie close to the positive real axis and negative imaginary axis
    LMP.sortquarticroots!(q)

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

    return SMatrix{2,2}(R11, R21, R12, R22), SMatrix{2,2}(dR11, dR21, dR12, dR22)
end
