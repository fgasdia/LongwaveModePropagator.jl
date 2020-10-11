#==
Functions related to identifying resonant modes ("eigenangles") within the
earth-ionosphere waveguide.

These functions use Budden's model, which represents the ground and ionosphere
as sharply reflecting boundaries. The middle of the waveguide is filled with a
fictitious medium with a refractive index that mimicks the propagation of the
radio wave through free space over curved earth.
==#

abstract type ModeEquation end

struct PhysicalModeEquation{W<:HomogeneousWaveguide,F} <: ModeEquation
    frequency::Frequency
    waveguide::W
    Mfcn::F
end

struct ModifiedModeEquation{W<:HomogeneousWaveguide,F} <: ModeEquation
    frequency::Frequency
    waveguide::W
    Mfcn::F
end

struct DZParams{F}
    ea::EigenAngle
    frequency::Frequency
    Mfcn::F
end

# TODO: G as its own type contained within this?
@with_kw struct SharpBoundaryVariables{T} @deftype T
    G12
    G32
    G33
    G11₁
    G13₁
    G31₁
    Δ₁
    Δ₁⁻¹
    P₁
    T₁
    G11₂
    G13₂
    G31₂
    Δ₂
    Δ₂⁻¹
    P₂
    T₂
    Δ
    Δ⁻¹
    q::MVector{4,T}
    B::SVector{5,T}
end

function Base.show(io::IO, ::MIME"text/plain", s::SharpBoundaryVariables{T}) where {T}
    println(io, "SharpBoundaryVariables{$T}: ")
    for n in fieldnames(SharpBoundaryVariables)
        if n != last(fieldnames(SharpBoundaryVariables))
            println(io, " $n: ", getfield(s, n))
        else
            # To avoid double newlines (one is automatic), use `print` on last field
            print(io, " $n: ", getfield(s, n))
        end
    end
end

@with_kw struct FresnelReflectionVariables{T<:Complex} @deftype T
    Ng²
    CNg²
    sqrtNg²mS²
    Rg::SDiagonal{2,T}
end

"""
    isroot(z::Complex; atol=1e-2)

Return `true` if both real and imaginary components of `z` are approximately equal to 0.
"""
function isroot(z::Complex; atol=1e-2)
    rz, iz = reim(z)
    isapprox(rz, 0, atol=atol) && isapprox(iz, 0, atol=atol)
end

"""
    isroot(x::Real; atol=1e-2)

Return `true` if the value of `x` is approximately equal to 0.
"""
isroot(x::Real; atol=1e-2) = isapprox(z, 0, atol=atol)

"""
    filterroots!(roots, frequency, waveguide; atol=1e-2)

Remove elements from `roots` if they are not valid roots (determined by `isroot`) to the modal
equation.

The absolute tolerance `atol` is passed through to `isroot`.
"""
function filterroots!(roots, frequency, waveguide; atol=1e-2)
    Mfcn(alt) = susceptibility(alt, frequency, waveguide)
    modeequation = PhysicalModeEquation(frequency, waveguide, Mfcn)

    i = 1
    while i <= length(roots)
        f = solvemodalequation(EigenAngle(roots[i]), modeequation)
        isroot(f, atol=atol) ? (i += 1) : deleteat!(roots, i)
    end
    return nothing
end

##########
# Reflection coefficients
##########

function _sharpboundaryreflection(ea::EigenAngle, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    q, B = bookerquartic(ea, M)

    # We choose the 2 roots corresponding to upward travelling waves as being
    # those that lie close to the positive real axis and negative imaginary axis
    q = sortquarticroots!(q)

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

    G11₂ = M11p1 - q[2]^2
    G13₂ = M[1,3] + q₂S
    G31₂ = M[3,1] + q₂S

    Δ₂ = G11₂*G33 - G13₂*G31₂
    Δ₂⁻¹ = 1/Δ₂
    P₂ = (-G12G33 + G13₂*G32)*Δ₂⁻¹
    T₂ = q[2]*P₂ - S*(-G11₂*G32 + G12*G31₂)*Δ₂⁻¹

    Δ = (T₁*C + P₁)*(C + q[2]) - (T₂*C + P₂)*(C + q[1])
    Δ⁻¹ = 1/Δ

    return SharpBoundaryVariables(G12, G32, G33, G11₁, G13₁, G31₁, Δ₁, Δ₁⁻¹, P₁, T₁,
                                  G11₂, G13₂, G31₂, Δ₂, Δ₂⁻¹, P₂, T₂, Δ, Δ⁻¹, q, B)
end

"""
    sharpboundaryreflection(ea, M)

Return ionosphere reflection matrix `R` for a sharply bounded anisotropic ionosphere.

[^Sheddy1968a] introduces the matrix equation ``G E = 0`` in the coordinate system of Budden
where the dispersion matrix ``G = I + M + L``. Nontrivial solutions to the equation require
``det(G) = 0``, which may be written as a quartic in ``q = n cos(θ)``. Elements of the
dispersion matrix are then used to calculate the reflection matrix `R`.

The reflection coefficient matrix for the sharply bounded case is used as a starting solution
for integration of the reflection coefficient matrix through the ionosphere.

# References

[^Sheddy1968a]: C. H. Sheddy, “A General Analytic Solution for Reflection From a Sharply Bounded Anisotropic Ionosphere,” Radio Science, vol. 3, no. 8, pp. 792–795, Aug. 1968.
"""
function sharpboundaryreflection(ea::EigenAngle, M)
    C = ea.cosθ
    C2 = 2*C

    @unpack P₁, T₁, P₂, T₂, Δ⁻¹, q = _sharpboundaryreflection(ea, M)

    T₁C = T₁*C
    T₂C = T₂*C

    R11 = ((T₁C - P₁)*(C + q[2]) - (T₂C - P₂)*(C + q[1]))*Δ⁻¹  # ∥R∥
    R22 = ((T₁C + P₁)*(C - q[2]) - (T₂C + P₂)*(C - q[1]))*Δ⁻¹  # ⟂R⟂
    R12 = -C2*(T₁*P₂ - T₂*P₁)*Δ⁻¹  # ⟂R∥
    R21 = -C2*(q[1] - q[2])*Δ⁻¹  # ∥R⟂

    return SMatrix{2,2}(R11, R21, R12, R22)
end

# TODO: Autodiff with Zygote?
function sharpboundaryreflection(ea::EigenAngle, M, ::Dθ)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    @unpack G12, G32, G33, G11₁, G13₁, G31₁, Δ₁, Δ₁⁻¹, P₁, T₁,
        G11₂, G13₂, G31₂, Δ₂, Δ₂⁻¹, P₂, T₂, Δ, Δ⁻¹, q, B = _sharpboundaryreflection(ea, M)

    # Precompute some variables
    C2 = 2*C

    T₁C = T₁*C
    T₂C = T₂*C

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

"""
    wmatrix(ea::EigenAngle, T)

Compute submatrix elements of `W` for solving the differential equation of reflection matrix
`R` wrt `z`.

``W`` is also known as ``S`` in many texts.

Following Budden's [^Budden1955a] formalism for the reflection matrix of a plane wave
obliquely incident on the ionosphere, the wave below the ionosphere can be resolved into
upgoing and downgoing waves of elliptical polarization, each of whose components are
themselves resolved into a component with the electric field in the plane of propagation and
a component perpendicular to the plane of propagation. The total field can be written in
matrix form as ``e = L q`` where ``L`` is a 4×4 matrix that simply selects and specifies the
incident angle of the components and ``q`` is a column matrix of the complex amplitudes of
the component waves. By inversion, ``q = L⁻¹ e`` and its derivative wrt height `z` is
``q′ = -i L⁻¹ T L q = -½ i S q``. Then ``S = 2 L⁻¹ T L``, which is calculated by this
function, describes the change in amplitude of the upgoing and downgoing component waves.

See also: [`tmatrix`](@ref), [`susceptibility`](@ref)

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
[^Sheddy1968]: C. H. Sheddy, R. Pappert, Y. Gough, and W. Moler, “A Fortran program for mode constants in an earth-ionosphere waveguide,” Naval Electronics Laboratory Center, San Diego, CA, Interim Report 683, May 1968.
"""
function wmatrix(ea::EigenAngle, T)
    C = ea.cosθ
    Cinv = ea.secθ

    # Precompute
    T12Cinv = T[1,2]*Cinv
    T14Cinv = T[1,4]*Cinv
    T32Cinv = T[3,2]*Cinv
    T34Cinv = T[3,4]*Cinv
    CT41 = C*T[4,1]

    #==
    Sheddy 1968 has at least 1 "typo": S11b should be written T14/C+CT41
    Budden 1988 has at least 1 typo: -T11 - T14 + CT41 + T44 for element [1,1] of [2,1]

    The analytical solution below was checked using SymPy.
    There is also a test against 2*inv(L)*T*L:

    T = SMatrix{4,4}(T11, 0, T31, T41,
                     T12, 0, T32, T42,
                     0, 1, 0, 0,
                     T14, 0, T34, T44)
    L = SMatrix{4,4}(C, 0, 0, 1,
                     0, -1, -C, 0,
                     -C, 0, 0, 1,
                     0, -1, C, 0)
    # 2*inv(L) = Linv2
    L2inv = SMatrix{4,4}(Cinv, 0, -Cinv, 0,
                        0, -1, 0, -1,
                        0, -Cinv, 0, Cinv,
                        1, 0, 1, 0)
    W = Linv2*T*L

    # W = 2*(L\T)*L

    ---

    W = | W11 | W12 |
        | W21 | W22 |

    W11 = | a11+a11r | -b11 |
          | -c11     | d11  |

    W12 = | a21+a21r | -b11 |
          | c12      | d12  |

    W21 = | a21-a21r | b21  |
          | c11      | -d12 |

    W22 = | a11-a11r | b21  |
          | -c12     | -d11 |
    ==#

    a11 = T[1,1] + T[4,4]
    a11r = T14Cinv + CT41
    a21 = T[4,4] - T[1,1]
    a21r = T14Cinv - CT41

    b11 = T12Cinv + T[4,2]
    b21 = T12Cinv - T[4,2]

    c11 = T[3,1] + T34Cinv
    c12 = T[3,1] - T34Cinv

    d11 = C + T32Cinv
    d12 = T32Cinv - C

    # Form the four 2x2 submatrices of `S`
    W11 = SMatrix{2,2}(a11+a11r, -c11, -b11, d11)
    W12 = SMatrix{2,2}(a21+a21r, c12, -b11, d12)
    W21 = SMatrix{2,2}(a21-a21r, c11, b21, -d12)
    W22 = SMatrix{2,2}(a11-a11r, -c12, b21, -d11)

    return W11, W21, W12, W22
end

"""
    dwmatrixdθ(ea, M, T)

Return the 4 submatrix elements of the derivative of the `W` matrix wrt θ.

See also: [`wmatrix`](@ref), [`susceptibility`](@ref), [`tmatrix`](@ref)
"""
function dwmatrixdθ(ea::EigenAngle, M, T)
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = ea.secθ
    C²inv = Cinv^2

    dS = C
    dC = -S
    dC² = -2*S*C
    dCinv = S*C²inv

    den = inv(1 + M[3,3])

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

    dt12Cinv = dT12*Cinv + T[1,2]*dCinv
    dt14Cinv = dT14*Cinv + T[1,4]*dCinv
    dt32Cinv = dT32*Cinv + T[3,2]*dCinv
    dt34Cinv = dT34*Cinv + T[3,4]*dCinv
    dt41C = dC*T[4,1]

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
    dW11 = SMatrix{2,2}(ds11a+ds11b, -ds21, -ds12, ds22)
    dW12 = SMatrix{2,2}(-dd11a+dd11b, dd21, -ds12, -dd22)
    dW21 = SMatrix{2,2}(-dd11a-dd11b, ds21, dd12, dd22)
    dW22 = SMatrix{2,2}(ds11a-ds11b, -dd21, dd12, -ds22)

    return dW11, dW21, dW12, dW22
end

"""
    dRdz(R, params, z)

Return the differential of the reflection matrix `R` wrt height `z`.

Following the Budden [^Budden1955a] formalism for the reflection of an (obliquely) incident
plane wave from a horizontally stratified ionosphere, the differential of the reflection
matrix `R` with height `z` can be described by ``2i R′ = W⁽²¹⁾ + W⁽²²⁾R - RW⁽¹¹⁾ - RW⁽¹²⁾R``.
Integrating ``R′`` wrt height `z`, gives the reflection matrix ``R`` for the ionosphere.

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
"""
function dRdz(R, params, z)
    @unpack ea, frequency, Mfcn = params

    k = frequency.k

    M = Mfcn(z)
    T = tmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)

    # the factor k/(2i) isn't explicitly in Budden1955a b/c of his change of variable ``s = kz``
    return -1im/2*k*(W21 + W22*R - R*W11 - R*W12*R)
end

function dRdθdz(RdRdθ, params, z)
    @unpack ea, frequency, Mfcn = params

    k = frequency.k

    M  = Mfcn(z)
    T = tmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)
    dW11, dW21, dW12, dW22 = dwmatrixdθ(ea, M, T)

    R = RdRdθ[SVector(1,2),:]
    dRdθ = RdRdθ[SVector(3,4),:]

    dz = -1im/2*k*(W21 + W22*R - R*W11 - R*W12*R)
    dθdz = -1im/2*k*(dW21 + dW22*R + W22*dRdθ - (dRdθ*W11 + R*dW11) -
                    (dRdθ*W12*R + R*dW12*R + R*W12*dRdθ))

    return vcat(dz, dθdz)
end

# TODO: Function args first
function integratedreflection(ea::EigenAngle, frequency::Frequency,
    waveguide::HomogeneousWaveguide, Mfcn::F) where {F}  # `F` forces dispatch on functions

    Mtop = susceptibility(TOPHEIGHT, frequency, waveguide)
    Rtop = sharpboundaryreflection(ea, Mtop)

    dzparams = DZParams(ea, frequency, Mfcn)
    prob = ODEProblem{false}(dRdz, Rtop, (TOPHEIGHT, BOTTOMHEIGHT), dzparams)

    #==
    | tolerance | method |
    |-----------|--------|
    | 1e-8      | Vern8  | slightly higher memory, much faster than Vern6, Vern7
    | 1e-6      | Vern8  |
    ==#

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, Vern8(), abstol=1e-8, reltol=1e-8,
                save_on=false, save_start=false, save_end=true)

    R = sol[end]

    return R
end

# This is kept as a completely separate function because the size of the matrix being
# integrated is different and therefore the size of sol[end] is different too
# The derivative terms are intertwined with the non-derivative terms so we can't do only
# the derivative terms
function integratedreflection(ea::EigenAngle, frequency::Frequency,
    waveguide::HomogeneousWaveguide, Mfcn::F, ::Dθ) where {F}

    Mtop = susceptibility(TOPHEIGHT, frequency, waveguide)
    RdRdθtop = sharpboundaryreflection(ea, Mtop, Dθ())

    dzparams = DZParams(ea, frequency, Mfcn)
    prob = ODEProblem{false}(dRdθdz, RdRdθtop, (TOPHEIGHT, BOTTOMHEIGHT), dzparams)

    # NOTE: When save_on=false, don't try interpolating the solution!
    # Purposefully higher tolerance than non-derivative version
    sol = solve(prob, Vern8(), abstol=1e-8, reltol=1e-8,
                save_on=false, save_start=false, save_end=true)

    R = sol[end]

    return R
end


##########
# Ground reflection coefficient matrix
##########

function _fresnelreflection(ea::EigenAngle, ground, frequency)
    C, S² = ea.cosθ, ea.sin²θ
    ω = frequency.ω

    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*E0))

    CNg² = C*Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²mS²)/(CNg² + sqrtNg²mS²)
    Rg22 = (C - sqrtNg²mS²)/(C + sqrtNg²mS²)

    Rg = SDiagonal(Rg11, Rg22)  # TODO: time this - slower than SMatrix?

    return FresnelReflectionVariables(Ng², CNg², sqrtNg²mS², Rg)
end

"""
    fresnelreflection(ea, ground, frequency)

Return the Fresnel reflection coefficient matrix for the ground free-space interface at the
ground (``z = 0``). Follows the formulation in [^Morfitt1976] pages 25-26.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency)
    @unpack Rg = _fresnelreflection(ea, ground, frequency)

    return Rg
end

function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency, ::Dθ)
    C, S, S² = ea.cosθ, ea.sinθ, ea.sin²θ
    ω = frequency.ω

    @unpack Rg, Ng², CNg², sqrtNg²mS² = _fresnelreflection(ea, ground, frequency)

    S2 = 2*S

    dRg11 = (S2*Ng²*(1 - Ng²))/(sqrtNg²mS²*(CNg² + sqrtNg²mS²)^2)
    dRg22 = (S2*(C - sqrtNg²mS²))/(sqrtNg²mS²*(sqrtNg²mS² + C))

    dRg = SDiagonal(dRg11, dRg22)

    return Rg, dRg
end

##########
# Identify EigenAngles
##########

"""
    modalequation(R, Rg)

Return value of the determinental mode equation ``det(Rg R - I)`` given `R` and `Rg`.

If the reflection matrix ``R₁`` represents the reflection coefficient at height ``z₁``, then
the reflection coefficient at ``z₂`` is ``R₂ = R₁ exp(2ik(z₂ - z₁) sin(θ))``. For a wave that
starts at the ground, reflects off the ionosphere at level ``h`` with ``R``, and then at the
ground with ``Rg``, the fields must be multiplied by ``̅Rg R exp(-2ikh sin(θ))``. For a
self-consistent waveguide mode, the resulting wave must be identifcal with the original
wave, and a necessary and sufficient condition for this is that one eigenvalue of this matrix
is unity. This leads to the mode equation [^Budden1962].

# References

[^Budden1962]: K. G. Budden and N. F. Mott, “The influence of the earth’s magnetic field on
radio propagation by wave-guide modes,” Proceedings of the Royal Society of London. Series A.
Mathematical and Physical Sciences, vol. 265, no. 1323, pp. 538–553, Feb. 1962.
"""
function modalequation(R, Rg)
    return det(Rg*R - I)
end

"""
See https://folk.ntnu.no/hanche/notes/diffdet/diffdet.pdf
"""
function modalequationdθ(R, dR, Rg, dRg)
    A = Rg*R - I
    dA = dRg*R + Rg*dR
    return det(A)*tr(A\dA)
end

function solvemodalequation(ea::EigenAngle, modeequation::PhysicalModeEquation)
    @unpack frequency, waveguide, Mfcn = modeequation
    R = integratedreflection(ea, frequency, waveguide, Mfcn)
    Rg = fresnelreflection(ea, waveguide.ground, frequency)

    f = modalequation(R, Rg)
    return f
end

"""
This returns R and Rg in addition to df because the only time this function is needed, we also
need R and Rg (in excitationfactors).
"""
function solvemodalequation(ea::EigenAngle, modeequation::PhysicalModeEquation, ::Dθ)
    # Dθ always uses `susceptibility` as Mfcn
    @unpack frequency, waveguide, Mfcn = modeequation

    RdR = integratedreflection(ea, frequency, waveguide, Mfcn, Dθ())
    R = RdR[SVector(1,2),:]
    dR = RdR[SVector(3,4),:]

    Rg, dRg = fresnelreflection(ea, waveguide.ground, frequency, Dθ())

    df = modalequationdθ(R, dR, Rg, dRg)
    return df, R, Rg
end

function findmodes(origcoords, grpfparams::GRPFParams, modeequation::ModeEquation)

    # WARNING: If tolerance of mode finder is much less than the R integration
    # tolerance, it may possible multiple identical modes will be identified. Checks for
    # valid and redundant modes help ensure valid eigenangles are returned from this function.

    @unpack frequency, waveguide = modeequation

    f(θ) = solvemodalequation(EigenAngle(θ), modeequation)
    roots, poles = grpf(f, origcoords, grpfparams)

    # Ensure roots are valid solutions to the modal equation
    filterroots!(roots, frequency, waveguide)

    # Remove any redundant modes
    # if tolerance is 1e-8, this rounds to 7 decimal places
    sort!(roots, by=reim)
    ndigits = round(Int, abs(log10(grpfparams.tolerance)+1), RoundDown)
    unique!(z->round(z, digits=ndigits), roots)

    return EigenAngle.(roots)
end

"""
    triangulardomain(Za, Zb, Zc, Δr)

Generate initial mesh node coordinates for a *grid-aligned right* triangle domain
∈ {`Za`, `Zb`, `Zc`} with initial mesh step `Δr`.

This function generates a mesh grid for particular right triangles where two
sides of the triangle are aligned to the underlying real/imaginary grid. Examples
of such triangles are:

a ----- b
|     /
|   /
| /
c

where
    - a, b have greatest extent in x
    - a, c have greatest extent in y
"""
function triangulardomain(Za::Complex, Zb::Complex, Zc::Complex, Δr)
    rZa, iZa = reim(Za)
    rZb, iZb = reim(Zb)
    rZc, iZc = reim(Zc)

    #==
    # Check if this is a right triangle
    validtriangle = true
    if rZa == rZb == rZc
        validtriangle = false
    elseif iZa == iZb == iZc
        validtriangle = false
    elseif rZa == rZb
        iZa == iZc || iZb == iZc || (validtriangle = false)
    elseif rZa == rZc
        iZa == iZb || iZb == iZc || (validtriangle = false)
    elseif rZb == rZc
        iZa == iZb || iZa == iZc || (validtriangle = false)
    else
        validtriangle = false
    end
    validtriangle || throw(ArgumentError("`Za`, `Zb`, `Zc` do not define a grid-aligned right triangle"))

    iZa == iZb || ((Zb, Zc) = (Zc, Zb))
    rZb > rZa || ((Za, Zb) = (Zb, Za))
    ==#

    # Determine `dx` and `dy`
    X = rZb - rZa
    Y = abs(iZa - iZc)

    n = ceil(Int, Y/Δr + 1)
    dy = Y/(n-1)
    half_dy = dy/2

    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    # precalculate
    mn = m*n
    slope = Y/X

    v = Vector{complex(typeof(dx))}()

    on = false
    for j = 0:m-1  # col
        for i = 0:n-1  # row (we're traversing down column)

            x = rZa + dx*j
            y = iZc + dy*i

            if (i+1) == n
                on = !on
            end
            if on
                y -= half_dy
            end

            if y >= (iZc + slope*(x - rZa))
                push!(v, complex(x, y))
            end
        end
    end

    return v
end

"""
    eiwgdomain(Zb, Ze, d, Δr)

This generates a vector of complex coordinates as a rectangulardomain with Δr on the left
and a uppertriangularrectdomain with Δr/5 on the right. The transition occurs at `d` in real.
"""
function eiwgdomain(Zb, Ze, d, Δr)
    if real(Zb) > d
        tricoords = uppertriangularrectdomain(Zb, Ze, Δr/5)
    else
        tricoords = uppertriangularrectdomain(complex(d, imag(Zb)), Ze, Δr/5)
        rectcoords = rectangulardomain(Zb, complex(d, imag(Ze)), Δr)
    end

    origcoords = vcat(tricoords, rectcoords)
    unique!(origcoords)

    return origcoords
end

function uppertriangularrectdomain(Zb::Complex, Ze::Complex, Δr)
    rZb, iZb = reim(Zb)
    rZe, iZe = reim(Ze)

    # Determine `dx` and `dy`
    X = rZe - rZb
    Y = abs(iZe - iZb)

    n = trunc(Int, cld(Y, Δr) + 1)
    dy = Y/(n-1)
    half_dy = dy/2

    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    slope = 1  # 45° angle

    v = Vector{complex(typeof(X))}()

    on = false
    for j = 0:m-1  # col
        x = rZb + dx*j

        for i = 0:n-1  # row (we're traversing down column)
            y = iZb + dy*i

            if (i+1) == n
                on = !on
            end
            if on
                y -= half_dy
            end

            if x <= (rZe - (iZe - y)/slope)
                push!(v, complex(x, y))
            end
        end
    end

    return v
end

# TODO: Finish this function
function targetdomain(modes::Vector{T}, delta::T, Δr) where {T<:Complex}
    v = Vector{T}(undef)
    for i in eachindex(modes)
        zb = v[i] - delta
        ze = v[i] + delta
        append!(rectangulardomain(zb, ze, Δr))
    end
    return v
end
