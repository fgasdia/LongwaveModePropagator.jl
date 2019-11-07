using SpecialFunctions
using LinearAlgebra
using StaticArrays
using DiffEqBase, OrdinaryDiffEq
using Parameters

using PolynomialRoots: roots!

using ModifiedHankelFunctionsOfOrderOneThird
using GRPF

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

abstract type FieldComponent end
struct Ez <: FieldComponent end
struct Ey <: FieldComponent end
struct Ex <: FieldComponent end

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

# TODO: Make this inherit from AbstractSparseArray
@with_kw struct TMatrix{T<:Number} @deftype T
    T11
    T12
    T14
    T31
    T32
    T34
    T41
    T42
    T44
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

    earthcurvature = 2(z₀ - z)/earthradius

    # Temporary variables
    U² = U^2
    Y² = Y^2

    YU = Y*U
    nYU = n*YU
    mYU = m*YU
    lYU = l*YU

    nY² = n*Y²
    lmY² = l*m*Y²
    lnY² = l*nY²
    mnY² = m*nY²

    # Elements of `M`
    M11 = U² - l²*Y²
    M21 = im*nYU - lmY²
    M31 = -im*mYU - lnY²
    M12 = -im*nYU - lmY²
    M22 = U² - m²*Y²
    M32 = im*lYU - mnY²
    M13 = im*mYU - lnY²
    M23 = -im*lYU - mnY²
    M33 =  U² - n²*Y²

    M = SMatrix{3,3}(M11, M21, M31,
                     M12, M22, M32,
                     M13, M23, M33)

    M *= -X/(U*(U² - Y²))

    M -= earthcurvature*I  # This correction only occurs once after adding species?

    return M
end

function tmatrix(ea::EigenAngle, M::AbstractArray)
    S, C² = ea.sinθ, ea.cos²θ

    den = 1/(1 + M[3,3])

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

    return TMatrix(T11, T12, T14, T31, T32, T34, T41, T42, T44)
end


"""
Compute matrix elements for solving differential of reflection matrix `R` wrt `z`.

See Budden 1955 second method.

```math
e′ = -iTe
```
"""
function smatrix(ea::EigenAngle, T::TMatrix)
    C = ea.cosθ
    Cinv = 1/C

    # Temporary matrix elements T
    @unpack T11, T12, T14, T31, T32, T34, T41, T42, T44 = T

    # based on Sheddy et al., 1968, A Fortran Program...
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

function dsmatrixdθ(ea::EigenAngle, M, T::TMatrix)
    # Unpack
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = 1/C
    C²inv = 1/C²

    dS = C
    dC = -S
    dC² = -2*S*C
    dCinv = S/C²

    # Temporary matrix elements T
    @unpack T12, T14, T32, T34, T41 = T

    den = 1/(1 + M[3,3])

    dT11 = -dS*M[3,1]*den
    dT12 = dS*M[3,2]*den
    # dT13 = 0
    dT14 = dC²*den
    # dT21 = 0
    # dT22 = 0
    # dT23 = 0
    # dT24 = 0
    # dT31 = 0
    dT32 = dC²
    # dT33 = 0
    dT34 = dS*M[2,3]*den
    # dT41 = 0
    # dT42 = 0
    # dT43 = 0
    dT44 = -dS*M[1,3]*den

    dt12Cinv = dT12*Cinv + T12*dCinv
    dt14Cinv = dT14*Cinv + T14*dCinv
    dt32Cinv = dT32*Cinv + T32*dCinv
    dt34Cinv = dT34*Cinv + T34*dCinv
    dt41C = dC*T41

    ds11a = dT11 + dT44
    dd11a = dT11 - dT44
    ds11b = dt14Cinv + dt41C
    dd11b = dt14Cinv - dt41C
    ds12 = dt12Cinv
    dd12 = dt12Cinv
    ds21 = dt34Cinv
    dd21 = -dt34Cinv
    ds22 = dC + dt32Cinv
    dd22 = dC - dt32Cinv

    # Form the four 2x2 submatrices of `dS`
    dS11 = @SMatrix [ds11a+ds11b -ds12;
                     -ds21 ds22]
    dS12 = @SMatrix [-dd11a+dd11b -ds12;
                     dd21 -dd22]
    dS21 = @SMatrix [-dd11a-dd11b dd12;
                     ds21 dd22]
    dS22 = @SMatrix [ds11a-ds11b dd12;
                     -dd21 -ds22]

    return dS11, dS12, dS21, dS22
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
    angq = rad2deg(angle(qval))
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
        S*(-dD11_1*D32 + D12*dD31_1)*invΔ_1 -
        S*(-D11_1*D32 + D12*D31_1)*dinvΔ_1

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

    # dR
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
    T = tmatrix(ea, M)
    S = smatrix(ea, T)

    return -im/2*k*(S[3] + S[4]*R - R*S[1] - R*S[2]*R)
end

function dRdθdz(RdRdθ, params, z)
    ω, k = params.ω, params.k
    ea = params.ea
    z0, species, bfield = params.referenceheight, params.species, params.bfield

    M = susceptibility(ω, z0, z, species, bfield)
    T = tmatrix(ea, M)
    S = smatrix(ea, T)
    dS = dsmatrixdθ(ea, M, T)

    R = RdRdθ[1:4]
    dRdθ = RdRθ[5:end]

    dz = -im/2*k*(S[3] + S[4]*R - R*S[1] - R*S[2]*R)
    dzdθ = -im/2*k*(dS[3] + dS[4]*R + S[4]*dR -
        (dRdθ*S[1] + R*dS[1]) -
        (dRdθ*S[2]*R + R*dS[2]*R + R*S[2]*dRdθ))

    return vcat(dz, dzdθ)
end

function integratethroughionosphere(
    ea::EigenAngle,
    source::AbstractSource,
    fromheight,
    toheight,
    referenceheight,
    species,
    bfield::BField;
    deriv=false
)
    ω, k, λ = source.ω, source.k, source.λ  # k in km

    M = susceptibility(ω, referenceheight, fromheight, species, bfield)

    params = (ω=ω, k=k, ea=ea, referenceheight=referenceheight, species=species, bfield=bfield)

    if deriv
        R0, dR0 = sharplybounded_R_dRdθ(ea, M)
        prob = ODEProblem{false}(dRdθdz, vcat(R0, dR0), (fromheight, toheight), params)
        sol = solve(prob, Tsit5())
    else
        R0 = sharplyboundedR(ea, M)

        prob = ODEProblem{false}(dRdz, R0, (fromheight, toheight), params)
        sol = solve(prob, Tsit5()) #, reltol=1e-6)#, dtmax=λ/20)
    end
end

"""
Fresnel reflection coefficients for the ground free-space interface at the ground (z=0).

From Morfitt Shellman 1976 pg 25 (eq 71 & 72)
"""
function fresnelreflection(ea::EigenAngle, source::AbstractSource, ground::Ground)
    C, S² = ea.cosθ, ea.sin²θ

    ng² = ground.ϵᵣ - im*ground.σ/(source.ω*ϵ₀)

    tmp1 = C*ng²
    tmp2 = sqrt(ng² - S²)

    Rg11 = (tmp1 - tmp2)/(tmp1 + tmp2)
    Rg22 = (C - tmp2)/(C + tmp2)

    # TODO: Custom type
    return SMatrix{2,2}(Rg11, 0, 0, Rg22)
end

function fresnelreflectiondθ(ea::EigenAngle, source::AbstractSource, ground::Ground)
    C, S, S² = ea.cosθ, ea.sinθ, ea.sin²θ

    ng² = ground.ϵᵣ - im*ground.σ/(source.ω*ϵ₀)

    tmp1 = C*ng²
    tmp2 = sqrt(ng² - S²)

    Rg11 = (tmp1 - tmp2)/(tmp1 + tmp2)
    Rg22 = (C - tmp2)/(C + tmp2)

    dRg11 = (2*ng²*(1 - ng²)*S)/(tmp2*(tmp1 + tmp2)^2)
    dRg22 = (2*(C - tmp2)*S)/(tmp2*(tmp2 + C))

    # TODO: Custom type
    return SMatrix{2,2}(Rg11, 0, 0, Rg22), SMatrix{2,2}(dRg11, 0, 0, dRg22)
end

"""
Determinental mode equation assuming `R` and `Rg` at θ

XXX: Maybe explicitly make this F(θ) and have it calculate Rg and R?
XXX: Maybe hardcode this?
"""
modalequation(R, Rg) = det(Rg*R .- 1)

"""
See https://folk.ntnu.no/hanche/notes/diffdet/diffdet.pdf
"""
function modalequationdθ(R, dR, Rg, dRg)
    A = Rg*R .- 1
    dA = dRg*R + Rg*dR
    return det(A)*tr(inv(A)*dA)
end

"""
See Morfitt 1980 Simplified ... Eqs 26-32

Morfitt specifies that `S` is the sine of θ at reference height `h` which is apparently(?)
not the same as his `b`, which is Wait's effective reflection height. But the `R`s and `Rg`s
are elements of the reflection matrix looking into iono and towards ground from the level `d`.
Therefore `T`s are also at `d`. BUT (on page 20) "the mode conversion program, as they are
now programmed, require that the height variable `d` be at the ground so that in equations
(26) through (32), `z = d = 0`".

TODO: These names are confusing - this function also returns height gains?

sw_wvgd.for
"""
function modeparameters(ea::EigenAngle)
    S = ea.sinθ

    sqrtS = sqrt(S)

    F₁, F₂, F₃, F₄ = heightgain()
    h₁qd, h₂qd, h₁′qd, h₂′qd = modhankel(qd)

    D11 = (F₁*h₁qd + F₂*h₂qd)^2
    D12 = (F₁*h₁qd + F₂*h₂qd)*(F₃*h₁qd + F₄h₂qd)
    D22 = (F₃h₁qd + F₄*h₂qd)^2

    T₁ = (sqrtS*(1 + Rg[1,1])^2*(1 - R[2,2]*Rg[2,2])) / (dFdθ*Rg[1,1]*D11)
    T₂ = (sqrtS*(1 + Rg[2,2])^2*(1 - R[1,1]*Rg[1,1])) / (dFdθ*Rg[2,2]*D22)
    T₃ = (sqrtS*(1 + Rg[1,1])*(1 + Rg[2,2])*R[2,1]) / (dFdθ*D12)
    T₄ = R[1,2]/R[2,1]

    τ₁ = D11*T₁
    τ₂ = D22*T₂
    τ₃ = D12*T₃

    # `f` = `Ey/Hy` at the ground
    # In LWPC (not sure where this requirement comes from). M 1980 says either (eq 33)
    if abs2(1 - R[1,1]*Rg[1,1]) > abs2(1 - R[2,2]*Rg[2,2])
        f = T2/(T3*T4)
    else
        f = T3/T1
    end

    # height gain for vertical electric field Ez
    f₁(z) = exp((z-d)/a)*(F₁*h₁q + F₂*h₂q)/(F₁*h₁q₀ + F₂*h₂q₀)
    # horizontal electric field Ex
    f₂(z) = 1/(im*k)*df₁dz  # TODO
    # Ey, normal to plane of propagation
    f₃(z) = (F₃*h₁q + F₄*h₂q)*f/(F₃*h₁q₀ + F₄*h₂q₀)

    if vertical
        λ_Ez = τ₁*S²
        λ_Ex = τ₁*S
        λ_Ey = -τ₃*S/f
    elseif endon
        λ_Ez = -τ₁*S
        λ_Ex = -τ₁
        λ_Ey = τ₃/f
    elseif broadside
        λ_Ez = -τ₃*T₄*S/f
        λ_Ex = -τ₃*T₄/f
        λ_Ey = τ₂/f^2
    end

    return f₁, f₂, f₃, λ_Ez, λ_Ex, λ_Ey
end

"""
Morfitt 1980

NOTE: LWPC and Morfitt 1980 have different conventions for (1,2,3). Morfitt has (1,2,3) → (z, x, y)
LWPC has (1,2,3) → (z, y, x) at least for heightgains...

E field strength amplitude in dB above μV/m for 1 kW radiator. Phase relative to free space

TODO: Special function for vertical Ez only, zt=zr=0, and γ=φ=0. (Like GRNDMC)
"""
function Efield()
    Q = 3.248e-5*sqrt(1000)*k/sqrt(f)  # Morfitt 1980 eq 41, adjusted for MKS units
    tmp = Q/sqrt(sin(x/a))

    # Transmit antenna orientation (electric dipole)
    # See Morfitt 1980 pg 22
    Sγ, Cγ = sincos(γ)
    Sφ, Cφ = sincos(φ)

    ρ = rcvrrange  # NOTE: don't actually multiply this on this function, so we can generalize
    # within the guide later

    modesum = zero()
    for ea in modes
        # calculate mode params. Note: All 3 directions needed for xmtr, only "wanted" term needed for rcvr
        xmtrterm = λ_Ez*f₁(zt)*Cγ + λ_Ex*f₂(zt)*Sγ*Cφ + λ_Ey*f₃(zt)*Sγ*Sφ
        rcvrterm = f₁(zr)  # choose which coordinate is wanted
        modsum += xmtrterm*rcvrterm*exp(-im*k)*(S - 1)*ρ
    end
end

"""
MS 76 (pg 38)
Pappert et al 1970 A FORTRAN program...
Pappert and Shockey 1971 WKB Mode Summing Program...

Morfitt 1980 Simplified Eqs 16-23
According to Morfitt, `H` (his `b`) is the effective ionospheric height. The value to be
assigned to this variable has been determined to be that height level for which the conductivity
parameter ωᵣ is ``2.5×10⁵`` sec⁻¹. This value of `b` is the reference height (`h'`) as
defined by Wait in terms of exponential profiles.

fvert = height gain for vertical electric field cponent (Ez)
fhorz = height gain for horizontal electric field component (Ey) normal to plane of propagation
g = height gain for horizontal electric field component (Ex) in plane of propagation
"""
function heightgain(z, ea::EigenAngle, ground::Ground)
    C², S² = ea.cos²θ, ea.sin²θ

    α = 2/earthradius
    tmp = (α/k)^(2/3)

    n₀² = 1 - α*(H - z)
    Ng² = ground.ϵᵣ - im*ground.σ/(ω*ϵ₀)

    q = h -> (C² - α*(H - h))/tmp
    q₀ = q(0)
    qd = q(d)

    h₁q₀, h₂q₀, h₁′q₀, h₂′q₀ = modhankel(q₀)
    h₁qd, h₂qd, h₁′qd, h₂′qd = modhankel(qd)

    H₁q₀ = h₁′q₀ + tmp*h₁q₀/2
    H₂q₀ = h₂′q₀ + tmp*h₂q₀/2

    H₁qd = h₁′qd + tmp*h₁qd/2
    H₂qd = h₂′qd + tmp*h₂qd/2

    t1 = n₀²/Ng²
    t2 = cbrt(k/α)
    t3 = sqrt(Ng²-S²)

    F₁ = -(H₂q₀ - im*t1*t2*t3*h₂q₀)
    F₂ = H₁q₀ - im*t1*t2*t3*h₁q₀
    F₃ = -(h₂′q₀ - im*t2*t3*h₂q₀)
    F₄ = h₁′q₀ - im*t2*t3*h₁q₀  # NOTE: Typo in P&S 71 eq 9

    return F₁, F₂, F₃, F₄

    # fvert = h -> exp((h-d)/earthradius)*(F₁*h₁(q(h)) + F₂*h₂(q(h)))/(F₁*h₁qd + F₂*h₂qd)
    # fhorz = h -> (F₃*h₁(q(h)) + F₄*h₂(q(h)))/(F₃*h₁qd + F₄*h₂qd)
    # g = h -> 1/(im*k)
end
