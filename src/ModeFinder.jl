using SpecialFunctions
using LinearAlgebra
using StaticArrays
using DiffEqBase, OrdinaryDiffEq
using Parameters

using PolynomialRoots: roots!

using ModifiedHankelFunctionsOfOrderOneThird

struct EigenAngle{T}
    θ::T  # radians, because df/dθ are in radians
    cosθ::T
    sinθ::T
    cos²θ::T
    sin²θ::T

    function EigenAngle{T}(θ::T) where T <: Number
        C = cos(θ)
        S = sin(θ)
        C² = C^2
        S² = 1 - C²
        new(θ, C, S, C², S²)
    end
end
EigenAngle(θ::T) where T <: Number = EigenAngle{T}(θ)

@enum FieldComponent begin
    Ez
    Ey
    Ex
end

@with_kw struct SharplyBoundedR{T<:Number} @deftype T
    q::MVector{4, ComplexF64}
    B::MVector{5, T}
    D12
    D32
    D33
    D11_1
    D13_1
    D31_1
    Δ_1
    invΔ_1
    P_1
    T_1
    D11_2
    D13_2
    D31_2
    Δ_2
    invΔ_2
    P_2
    T_2
    Δ
    invΔ
end

"""
Search boundary for zeros in complex plane.

TODO: This seems a little arbitrary...

See also: `lwp_input.for`
"""
function boundaries(freq)
    if freq < 20e3
        Zb = complex(60, 0)
        Ze = complex(90, -9)
        return Zb, Ze
    else
        Zb = complex(40, 0)
        Ze = complex(90, -6)
        return Zb, Ze
    end
end

"""
Computation of susceptibility `M` matrix as defined by Budden (1955)
[](10.1098/rspa.1955.0027).

`z` is "current" height
`z₀` is reference height for earth curvature where the index of refraction ≈ 1

# TODO: Add M up for each species
"""
function susceptibility(ω, z₀, z, spec::Constituent, bfield::BField)
    # Unpack
    B, l, m, n = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    l², m², n² = l^2, m^2, n^2

    # Constitutive relations (see Budden 1955, pg. 517)
    e, m, N, ν = spec.charge, spec.mass, spec.numberdensity, spec.collisionfrequency
    X = N(z)*e^2/(ϵ₀*m*ω^2)
    Y = e*B/(m*ω)
    Z = ν(z)/ω
    U = 1 - im*Z

    U² = U^2
    Y² = Y^2

    earthcurvature = 2(z₀ - z)/earthradius

    # In LWPC and Sheddy 1968 Fortran Program `earthcurvature` is not multiplied by capd
    # (the above line), even though it _is_ multiplied in MS 1976.
    # This seems to be supported by Pappert 1968
    M11 = U² - l²*Y²
    M21 = im*n*Y*U - l*m*Y²
    M31 = -im*m*Y*U - l*n*Y²
    M12 = -im*n*Y*U - l*m*Y²
    M22 = U² - m²*Y²
    M32 = im*l*Y*U - m*n*Y²
    M13 = im*m*Y*U - l*n*Y²
    M23 = -im*l*Y*U - m*n*Y²
    M33 =  U² - n²*Y²

    M = SMatrix{3,3}(M11, M21, M31,
                     M12, M22, M32,
                     M13, M23, M33)

    M *= -X/(U*(U² - Y²))

    M -= earthcurvature*I  # This correction only occurs once after adding species?

    return M
end

"""
Compute matrix elements for solving differential of reflection matrix `R` wrt `z`.

See Budden 1955 second method.

```math
e′ = -iTe
```
"""
function smatrix(ea::EigenAngle, M)
    # Unpack
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = 1/C

    den = 1/(1 + M[3,3])

    # Temporary matrix elements T
    m31d = M[3,1]*den
    m32d = M[3,2]*den

    T11 = -S*m31d
    T12 = S*m32d
    # T13 = 0
    T14 = (C² + M[3,3])*den
    # T21 = 0
    # T22 = 0
    # T23 = 1
    # T24 = 0
    T31 = M[2,3]*m31d - M[2,1]
    T32 = C² + M[2,2] - M[2,3]*m32d
    # T33 = 0
    T34 = S*M[2,3]*den
    T41 = 1 + M[1,1] - M[1,3]*m31d
    T42 = M[1,3]*m32d - M[1,2]
    # T43 = 0
    T44 = -S*M[1,3]*den

    # S matrix, based on Sheddy et al., 1968, A Fortran Program...
    t12Cinv = T12*Cinv
    t14Cinv = T14*Cinv
    t32Cinv = T32*Cinv
    t34Cinv = T34*Cinv
    t41C = C*T41

    s11a = T11 + T44
    d11a = T11 - T44
    s11b = t14Cinv + t41C
    d11b = t14Cinv - t41C
    s12 = t12Cinv + T42
    d12 = t12Cinv - T42
    s21 = T31 + t34Cinv
    d21 = T31 - t34Cinv
    s22 = C + t32Cinv
    d22 = C - t32Cinv

    # Form the four 2x2 submatrices of `S`
    S11 = @SMatrix [s11a+s11b -s12;
                    -s21 s22]
    S12 = @SMatrix [-d11a+d11b -s12;
                    d21 -d22]
    S21 = @SMatrix [-d11a-d11b d12;
                    s21 d22]
    S22 = @SMatrix [s11a-s11b d12;
                    -d21 -s22]

    return S11, S12, S21, S22

    #== Equivalent to (but faster than):
    T = SMatrix{4,4}(T11, T21, T31, T41,
                     T12, T22, T32, T42,
                     T13, T23, T33, T43,
                     T14, T24, T34, T44)

    L = SMatrix{4,4}(C, 0, 0, 1,
                     0, -1, -C, 0,
                     -C, 0, 0, 1,
                     0, -1, C, 0)

    Linv = SMatrix{4,4}(Cinv, 0, -Cinv, 0,
                        0, -1, 0, -1,
                        0, -Cinv, 0, Cinv,
                        1, 0, 1, 0)

    S = Linv*T*L
    ==#
end

"""
Calculation of Booker quartic for solution of `R` for a sharply bounded ionosphere.

Based on Sheddy 1968 A General Analytic Solution for Reflection from a Sharply...

See also: [`sharplyboundedX!`](@ref)
"""
function bookerquartic(ea::EigenAngle, M)
    # Initialize
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # Booker quartic coefficients
    B4 = 1 + M[3,3]
    B3 = S*(M[1,3] + M[3,1])
    B2 = -(C² + M[3,3])*(1 + M[1,1]) + M[1,3]*M[3,1] -
         (1 + M[3,3])*(C² + M[2,2]) + M[2,3]*M[3,2]
    B1 = S*(M[1,2]*M[2,3] + M[2,1]*M[3,2] -
         (C² + M[2,2])*(M[1,3] + M[3,1]))
    B0 = (1 + M[1,1])*(C² + M[2,2])*(C² + M[3,3]) +
         M[1,2]*M[2,3]*M[3,1] + M[1,3]*M[2,1]*M[3,2] -
         M[1,3]*(C² + M[2,2])*M[3,1] -
         (1 + M[1,1])*M[2,3]*M[3,2] -
         M[1,2]*M[2,1]*(C² + M[3,3])

    q = @MVector zeros(ComplexF64, 4)  # algorithm requires complex type
    B = MVector{5}(B0, B1, B2, B3, B4)

    roots!(q, B, NaN, 4, false)

    return q, B
end


"""
Calculate angle from 315°.

Used to sort in [`sharplyboundedX!`](@ref).
"""
function anglefrom315(qval)
    angq = rad2deg(angle(qval))  # rad2deg adds a few ns (negligible)
    angq < 0 && (angq += 360)
    angq < 135 && (angq += 360)
    return abs(angq - 315)
end

function _common_sharplyboundedR(ea::EigenAngle, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # Solve equation `D = 0` as a quartic
    q, B = bookerquartic(ea, M)

    # For high in the ionosphere, we choose 2 solutions that lie close to positive real and
    # negative imaginary axis (315° on the complex plane)
    sort!(q, by=anglefrom315)

    # Constant entries of dispersion matrix `D`
    D12 = M[1,2]
    D32 = M[3,2]
    D33 = C² + M[3,3]

    # Values for two solutions of the Booker quartic corresponding to upgoing waves
    D11_1 = 1 + M[1,1] - q[1]^2
    D13_1 = M[1,3] + q[1]*S
    D31_1 = M[3,1] + q[1]*S

    Δ_1 = D11_1*D33 - D13_1*D31_1
    invΔ_1 = 1/Δ_1
    P_1 = (-D12*D33 + D13_1*D32)*invΔ_1
    T_1 = q[1]*P_1 - S*(-D11_1*D32 + D12*D31_1)*invΔ_1

    D11_2 = 1 + M[1,1] - q[2]^2
    D13_2 = M[1,3] + q[2]*S
    D31_2 = M[3,1] + q[2]*S

    Δ_2 = D11_2*D33 - D13_2*D31_2
    invΔ_2 = 1/Δ_2
    P_2 = (-D12*D33 + D13_2*D32)*invΔ_2
    T_2 = q[2]*P_2 - S*(-D11_2*D32 + D12*D31_2)*invΔ_2

    # Computation of entries of reflection matrix `R`
    Δ = (T_1*C + P_1)*(C + q[2]) - (T_2*C + P_2)*(C + q[1])
    invΔ = 1/Δ

    return SharplyBoundedR(q, B, D12, D32, D33, D11_1, D13_1, D31_1, Δ_1, invΔ_1, P_1, T_1,
        D11_2, D13_2, D31_2, Δ_2, invΔ_2, P_2, T_2, Δ, invΔ)
end

"""
From Sheddy 1968, A General Analytic Solution for Reflection ...
"""
function sharplybounded_R(ea::EigenAngle, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    @unpack q, P_1, T_1, P_2, T_2, invΔ = _common_sharplyboundedR(ea, M)

    R11 = ((T_1*C - P_1)*(C + q[2]) - (T_2*C - P_2)*(C + q[1]))*invΔ  # ∥R∥
    R22 = ((T_1*C + P_1)*(C - q[2]) - (T_2*C + P_2)*(C - q[1]))*invΔ  # ⟂R⟂
    R12 = -2*C*(T_1*P_2 - T_2*P_1)*invΔ  # ⟂R∥
    R21 = -2*C*(q[1] - q[2])*invΔ  # ∥R⟂

    return SMatrix{2,2}(R11, R21, R12, R22)
end

function sharplybounded_R_dRdθ(ea::EigenAngle, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    @unpack q, B, D12, D32, D33, D11_1, D13_1, D31_1, Δ_1, invΔ_1, P_1, T_1,
        D11_2, D13_2, D31_2, Δ_2, invΔ_2, P_2, T_2, Δ, invΔ = _common_sharplyboundedR(ea, M)

    # Additional calculations required for dR/dθ
    dS = C
    dC = -S
    dC² = -2*S*C
    dB3 = dS*(M[1,3] + M[3,1])
    dB2 = -dC²*(2 + M[1,1] + M[3,3])
    dB1 = dS/S*B[2] - S*dC²*(M[1,3] + M[3,1])
    dB0 = dC²*(2*C²*(1 + M[1,1]) + M[3,3] + M[2,2] + M[1,1]*(M[3,3] + M[2,2]) -
        M[1,3]*M[3,1] - M[1,2]*M[2,1])

    dq_1 = -(((dB3*q[1] + dB2)*q[1] + dB1)*q[1] + dB0) /
        (((4*B[5]*q[1] + 3*B[4])*q[1] + 2*B[3])*q[1] + B[2])
    dq_2 = -(((dB3*q[2] + dB2)*q[2] + dB1)*q[2] + dB0) /
        (((4*B[5]*q[2] + 3*B[4])*q[2] + 2*B[3])*q[2] + B[2])

    dD33 = dC²

    dD11_1 = -2*q[1]*dq_1
    dD13_1 = dq_1*S + q[1]*dS
    dD31_1 = dD13_1  # dq_1*S + q[1]*dS

    dΔ_1 = dD11_1*D33 + D11_1*dD33 - dD13_1*D31_1 - D13_1*dD31_1
    dinvΔ_1 = -dΔ_1/Δ_1^2

    dP_1 = (-D12*dD33 + dD13_1*D32)*invΔ_1 + (D13_1*D32 - D12*D33)*dinvΔ_1
    dT_1 = dq_1*P_1 + q[1]*dP_1 -
        dS*(-D11_1*D32 + D12*D31_1)*invΔ_1 -
        S*((-dD11_1*D32 + D12*dD31_1)*invΔ_1 + (D11_1*D32 - D12*D31_1)*dinvΔ_1)

    dD11_2 = -2*q[2]*dq_2
    dD13_2 = dq_2*S + q[2]*dS
    dD31_2 = dD13_2  # dq_2*S + q[2]*dS

    dΔ_2 = dD11_2*D33 + D11_2*dD33 - dD13_2*D31_2 - D13_2*dD31_2
    dinvΔ_2 = -dΔ_2/Δ_2^2

    dP_2 = (-D12*dD33 + dD13_2*D32)*invΔ_2 + (D13_2*D32 - D12*D33)*dinvΔ_2
    dT_2 = dq_2*P_2 + q[2]*dP_2 -
        dS*(-D11_2*D32 + D12*D31_2)*invΔ_2 -
        S*(-dD11_2*D32 + D12*dD31_2)*invΔ_2 -
        S*(-D11_2*D32 + D12*D31_2)*dinvΔ_2

    dΔ = dT_1*C² + T_1*dC² + dT_1*C*q[2] + T_1*dC*q[2] + T_1*C*dq_2 + dP_1*C + P_1*dC +
        dP_1*q[2] + P_1*dq_2 - (dT_2*C² + T_2*dC²) -
        (dT_2*C*q[1] + T_2*dC*q[1] + T_2*C*dq_1) -
        (dP_2*C + P_2*dC) - (dP_2*q[1] + P_2*dq_1)
    dinvΔ = -dΔ/Δ^2

    # R
    R11 = ((T_1*C - P_1)*(C + q[2]) - (T_2*C - P_2)*(C + q[1]))*invΔ  # ∥R∥
    R22 = ((T_1*C + P_1)*(C - q[2]) - (T_2*C + P_2)*(C - q[1]))*invΔ  # ⟂R⟂
    R12 = -2*C*(T_1*P_2 - T_2*P_1)*invΔ  # ⟂R∥
    R21 = -2*C*(q[1] - q[2])*invΔ  # ∥R⟂

    # dR11 = invΔ*(dT_1*C² + T_1*dC² + dT_1*C*q[2] + T_1*dC*q[2] + T_1*dC*q[2] + T_1*C*dq_2 -
    #     (dP_1*C + P_1*dC) - (dP_1*q[2] + P_1*dq_2) - (dT_2*C² + T_2*dC²) -
    #     (dT_2*C*q[1] + T_2*dC*q[1] + T_2*C*dq_1) +
    #     (dP_2*C + P_2*dC) + (dP_2*q[1] + P_2*dq_1)) +
    #     dinvΔ*(T_1*C² + T_1*C*q[2] - P_1*C - P_1*q[2] - T_2*C² - T_2*C*q[1] + P_2*C + P_2*q[1])

    # dR11 = dinvΔ*(C²*(T_1 - T_2) + C*(T_1*q[2] - T_2*q[1] + P_2 - P_1) + P_2*q[1] - P_1*q[2]) +
    #     invΔ*(C²*(dT_1 - dT_2) + dC²*(T_1 - T_2) + dC*(T_1*q[2] - P_1 - T_2*q[1] + P_2) +
    #         C*(dT_1*q[2] + T_1*dq_2 + dP_2 - dP_1 - dT_2*q[1] - T_2*dq_1) +
    #         dP_2*q[1] + P_2*dq_1 - dP_1*q[2] - P_1*dq_2)
    dR11 = dinvΔ*R11*Δ +
        invΔ*(C²*(dT_1 - dT_2) + dC²*(T_1 - T_2) + dC*(T_1*q[2] - P_1 - T_2*q[1] + P_2) +
            C*(dT_1*q[2] + T_1*dq_2 + dP_2 - dP_1 - dT_2*q[1] - T_2*dq_1) +
            dP_2*q[1] + P_2*dq_1 - dP_1*q[2] - P_1*dq_2)
    dR12 = dinvΔ*R12*Δ + dC/C*R12 - 2*C*(dT_1*P_2 + T_1*dP_2 - dT_2*P_1 - T_2*dP_1)*invΔ
    dR21 = dinvΔ*R21*Δ + dC/C*R21 - 2*C*(dq_1 - dq_2)*invΔ
    dR22 = dinvΔ*R22*Δ +
        invΔ*(C²*(dT_1 - dT_2) + dC²*(T_1 - T_2) + dC*(T_2*q[1] - T_1*q[2] + P_1 - P_2) +
            C*(dP_1 - dT_1*q[2] + dT_2*q[1] - T_1*dq_2 + T_2*dq_1 - dP_2) +
            dP_2*q[1] + P_2*dq_1 - dP_1*q[2] - P_1*dq_2)

    return SMatrix{2,2}(R11, R21, R12, R22), SMatrix{2,2}(dR11, dR21, dR12, dR22)
end

"""
Calculate the derivative of the reflection matrix `R` wrt height `z`.

Expects S to be a tuple of `(S11, S12, S21, S22)`.

The factor of `k` appears explicitly because Budden 1955 derives R′ wrt a height
variable ``s'' which includes k.
"""
function dRdz(R, params, z)
    # TODO: BUG: This function needs a redesign, maybe don't use `params`
    # The function takes more than twice as long as the sum of its components

    ω, k = params.ω, params.k
    ea = params.ea
    z0, species, bfield = params.referenceheight, params.species, params.bfield

    M = susceptibility(ω, z0, z, species, bfield)
    S = smatrix(ea, M)

    return -im/2*k*(S[3] + S[4]*R - R*S[1] - R*S[2]*R)
end

function integratethroughionosphere(
    ea::EigenAngle,
    source::AbstractSource,
    fromheight,
    toheight,
    referenceheight,
    species,
    bfield::BField
)
    ω, k, λ = source.ω, source.k, source.λ  # k in km

    M = susceptibility(ω, referenceheight, fromheight, species, bfield)
    R0 = sharplyboundedR(ea, M)

    params = (ω=ω, k=k, ea=ea, referenceheight=referenceheight, species=species, bfield=bfield)
    prob = ODEProblem{false}(dRdz, R0, (fromheight, toheight), params)
    sol = solve(prob, Tsit5(), dtmax=λ/50_000) #, reltol=1e-6)#, dtmax=λ/20)
end

"""
Fresnel reflection coefficients for the ground free-space interface at the ground (z=0).

From Morfitt Shellman 1976 pg 25 (eq 71)
"""
function fresnelreflection(ea::EigenAngle, source::AbstractSource, ground::Ground)
    C, S² = ea.cosθ, ea.sin²θ

    ng² = ground.ϵᵣ - im*ground.σ/(source.ω*ϵ₀)

    tmp1 = C*ng²
    tmp2 = sqrt(ng² - S²)

    Rg11 = (tmp1 - tmp2)/(tmp1 + tmp2)
    Rg22 = (C - tmp2)/(C + tmp2)

    return SMatrix{2,2}(Rg11, 0, 0, Rg22)
end

"""
Determinental mode equation assuming `R` and `Rg` at θ

XXX: Maybe explicitly make this F(θ) and have it calculate Rg and R?
XXX: Maybe hardcode this?
"""
F(R, Rg) = det(Rg*R .- 1)

"""
Adjugate of A

https://github.com/JuliaDiff/ForwardDiff.jl/issues/197

# TODO: Hardcode for 2x2 or otherwise a more efficient implementation?
"""
adj(A) = det(A)*inv(A)

"""
Height gain terms from MS 76 (pg 38)

at height `z` with reference height `H`. `d` is reference height for solving F(θ) = 0
"""
function heightgains(z, ea::EigenAngle)
    C² = ea.cos²θ

    # TODO: Rather than calling these functions, just insert `ζd` into equations for fpar, fperp, g

    ζ = (k/α)^(2/3)*(C² + α*(z - H))
    ζd = (k/α)^(2/3)*(C² + α*(d - H))

    # height gain for vertical electric field cponent (Ez)
    fpar(z) = (QC1*h1(ζ) + GC*h2(ζ))*exp((z - d)/earthradius)

    # height gain for horizontal electric field component (Ey) normal to plane of propagation
    fperp(z) = AC2*h1(ζ) + BC2*h2(ζ)

    # height gain for horizontal electric field component (Ex) in plane of propagation
    #  ``1/(im*k) * d/dz(f∥(z))
    g(z) = -im*exp((z-d)/a) * (cbrt(2/(a*k))*(QC1*h1′(ζ) + Gc1*h2′(ζ)) +
        (2/(a*k))*(QC1*h1(ζ) + GC1*h2(ζ)))

    # We calculate D11, D12, and D22 because pairs of these are needed in excitation factor
    D11 = fpar(ζd)^2
    D12 = fpar(ζd)*fperp(ζd)
    D22 = fperp(ζd)^2

    return D11, D12, D22
end

"""
Excitation factors from MS 76 (pg 37)
"""
function excitationfactor(ea::EigenAngle, R, Rg, component::FieldComponent, source::Source)
    θ, S = ea.θ, ea.sinθ

    B1 = S^(5/2)/dFdθ(θ)
    B2 = -B1/S

    D11, D12, D22 = heightgains(z)  # z should be `d` which is height we evaluate F(θ) at

    if component == Ez
        if source.exciter == vertical
            return B1*(1 + Rg11)^2*(1 - Rg22*R22)/(Rg11*D11)
        elseif source.exciter == horizontal_endon
            return B2*(1 + Rg11)^2*(1 - Rg22*R22)/(Rg11*D11)
        elseif source.exciter == horizontal_broadside
            return B2*R11*(1 + Rg22)*(1 + Rg11)/D12
        else
            error("Source exciter not supported")
        end
    elseif component == Ey
        if source.exciter == vertical
            return -B1/S*R21*(1 + Rg11)*(1 + Rg22)/D12
        elseif source.exciter == horizontal_endon
            return -B2/S*R21*(1 + Rg11)*(1 + Rg22)/D12
        elseif source.exciter == horizontal_broadside
            return -B2/S*(1 + Rg22)^2*(1 - Rg11*R11)/(Rg22*D22)
        else
            error("Source exciter not supported")
        end
    elseif component == Ex
        if source.exciter == vertical
            return B1/S*(1 + Rg11)^2*(1 - Rg22*R22)/(Rg11*D11)
        elseif source.exciter == horizontal_endon
            return B2/S*(1 + Rg11)^2*(1 - Rg22*R22)/(Rg11*D11)
        elseif source.exciter == horizontal_broadside
            return B2/S*R12*(1 + Rg22)*(1 + Rg11)/D12
        else
            error("Source exciter not supported")
        end
    else
        error("FieldComponent not supported")
    end
end
