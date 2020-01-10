using SpecialFunctions
using LinearAlgebra
using StaticArrays
using DiffEqBase, OrdinaryDiffEq
using Parameters
using GRPF

using PolynomialRoots: roots!

using ModifiedHankelFunctionsOfOrderOneThird

const BOTTOMHEIGHT = 0.0  # should this be an actual const? Nothing works if it's not 0...
const TOPHEIGHT = 91e3  # temporary - should be part of an actual IntegrationParameters

struct Derivative end

@with_kw struct WaveguideParameters
    bfield::BField
    ground::Ground
    species::Constituent
end

@with_kw struct HeightGainConstants{T<:Number} @deftype T
    h₁q₀
    h₂q₀
    h₁′q₀
    h₂′q₀
    F₁
    F₂
    F₃
    F₄
end

@with_kw struct SharpBoundary{T<:Number} @deftype T
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

##########
# Reflection coefficient matrix from a sharply bounded ionosphere
##########

"""
    bookerquartic(ea, M)

Return roots `q` and the coefficients `B` of the Booker quartic.

The Booker quartic is used in the solution of `R` for a sharply bounded ionosphere. This
function uses the `PolynomialRoots` package to find the roots.

See also: [`sharplybounded_R`](@ref)

# References

[^Sheddy1968a]: C. H. Sheddy, “A General Analytic Solution for Reflection From a Sharply Bounded Anisotropic Ionosphere,” Radio Science, vol. 3, no. 8, pp. 792–795, Aug. 1968.
"""
function bookerquartic(ea::EigenAngle{T}, M) where T<:Number
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

    q = @MVector zeros(complex(T), 4)  # algorithm requires complex type
    B = MVector{5}(B0, B1, B2, B3, B4)

    roots!(q, B, NaN, 4, false)

    return q, B
end

"""
    anglefrom315(v)

Calculate the angle of `v` in radians from 315° on the complex plane.
"""
function anglefrom315(v::Complex)
    ang = rad2deg(angle(v))
    ang < 0 && (ang += 360)
    ang < 135 && (ang += 360)
    return abs(ang - 315)
end

function anglefrom315(v::Real)
    return oftype(v, 45)
end

function _common_sharplyboundedreflection(ea::EigenAngle{<:Number}, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # Solve equation `D = 0` as a quartic
    q, B = bookerquartic(ea, M)

    # For high in the ionosphere, we choose 2 solutions that lie close to positive real and
    # negative imaginary axis (315° on the complex plane)
    sort!(q, by=anglefrom315)

    # Constant entries of dispersion matrix `D`
    q1S = q[1]*S
    q2S = q[2]*S

    M11p1 = 1 + M[1,1]

    D12 = M[1,2]
    D32 = M[3,2]
    D33 = C² + M[3,3]

    # Values for two solutions of the Booker quartic corresponding to upgoing waves
    D11_1 = M11p1 - q[1]^2
    D13_1 = M[1,3] + q1S
    D31_1 = M[3,1] + q1S

    Δ_1 = D11_1*D33 - D13_1*D31_1
    invΔ_1 = inv(Δ_1)
    P_1 = (-D12*D33 + D13_1*D32)*invΔ_1
    T_1 = q[1]*P_1 - S*(-D11_1*D32 + D12*D31_1)*invΔ_1

    D11_2 = M11p1 - q[2]^2
    D13_2 = M[1,3] + q2S
    D31_2 = M[3,1] + q2S

    Δ_2 = D11_2*D33 - D13_2*D31_2
    invΔ_2 = inv(Δ_2)
    P_2 = (-D12*D33 + D13_2*D32)*invΔ_2
    T_2 = q[2]*P_2 - S*(-D11_2*D32 + D12*D31_2)*invΔ_2

    # Computation of entries of reflection matrix `R`
    Δ = (T_1*C + P_1)*(C + q[2]) - (T_2*C + P_2)*(C + q[1])
    invΔ = inv(Δ)

    return SharpBoundary(q, B, D12, D32, D33, D11_1, D13_1, D31_1, Δ_1, invΔ_1, P_1, T_1,
        D11_2, D13_2, D31_2, Δ_2, invΔ_2, P_2, T_2, Δ, invΔ)
end

"""
    sharplyboundedreflection(ea, M)

Return ionosphere reflection matrix `R` for a sharply bounded anisotropic ionosphere.

[^Sheddy1968a] introduces the matrix equation ``D E = 0`` in the coordinate system of Budden
where the dispersion matrix ``D = I + M + L``. Nontrivial solutions to the equation require
``det(D) = 0``, which may be written as a quartic in ``q = n cos(θ)``. Elements of the
dispersion matrix are then used to calculate the reflection matrix `R`.

The reflection coefficient matrix for the sharply bounded case is used as a starting solution
for integration of the reflection coefficient matrix through the ionosphere.

# References

[^Sheddy1968a]: C. H. Sheddy, “A General Analytic Solution for Reflection From a Sharply Bounded Anisotropic Ionosphere,” Radio Science, vol. 3, no. 8, pp. 792–795, Aug. 1968.
"""
function sharplyboundedreflection(ea::EigenAngle{<:Number}, M)
    C = ea.cosθ
    C2 = 2*C

    @unpack q, P_1, T_1, P_2, T_2, invΔ = _common_sharplyboundedreflection(ea, M)

    R11 = ((T_1*C - P_1)*(C + q[2]) - (T_2*C - P_2)*(C + q[1]))*invΔ  # ∥R∥
    R22 = ((T_1*C + P_1)*(C - q[2]) - (T_2*C + P_2)*(C - q[1]))*invΔ  # ⟂R⟂
    R12 = -C2*(T_1*P_2 - T_2*P_1)*invΔ  # ⟂R∥
    R21 = -C2*(q[1] - q[2])*invΔ  # ∥R⟂

    return SMatrix{2,2}(R11, R21, R12, R22)
end

function sharplyboundedreflection(ea::EigenAngle{<:Number}, M, ::Derivative)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    @unpack q, B, D12, D32, D33, D11_1, D13_1, D31_1, Δ_1, invΔ_1, P_1, T_1,
        D11_2, D13_2, D31_2, Δ_2, invΔ_2, P_2, T_2, Δ, invΔ = _common_sharplyboundedreflection(ea, M)

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

    # return SMatrix{2,2}(R11, R21, R12, R22), SMatrix{2,2}(dR11, dR21, dR12, dR22)
    return @SMatrix [R11 R12;
                     R21 R22;
                     dR11 dR12;
                     dR21 dR22]
end

##########
# Reflection coefficient matrix for a vertically stratified ionosphere
##########

"""
    susceptibility(z, ω, bfield, species)

Computation of susceptibility matrix `M` as defined by [^Budden1955a] method 2.

TODO: Handle multiple species

TODO: additional references on curvature correction
Includes first order correction for earth curvature by means of a fictitious refractive index.

Budden's formalism [^Budden1955a] assumes that the ionosphere is horizontally stratified and contains
electrons whose motion is influenced by the earth's magnetic field and is damped by collisions.
A plane wave of specified polarization is incident at some (generally oblique) angle from
below and we wish to find the amplitude, phase, and polarization of the reflected wave.

The differential equations for calculating the reflection coefficient are derived from
Maxwell's equations together with the constitutive relations for the ionosphere. The
constitutive relations for the ionosphere may be written ``P = ϵ₀ M ⋅ E`` where ``P`` is the
electric polarization vector, ``E`` is the electric field vector, and ``M`` is a 3×3
susceptibility matrix calculated by this function.

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
"""
function susceptibility(z, frequency::Frequency, bfield::BField, species::Constituent)
    B, l, m, n = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    l², m², n² = l^2, m^2, n^2
    ω = frequency.ω

    # Constitutive relations (see Budden1955a, pg. 517)
    e, m, N, ν = species.charge, species.mass, species.numberdensity, species.collisionfrequency
    mω = m*ω
    X = N(z)*e^2/(ϵ₀*mω*ω)  # ωₚ²/ω² plasma frequency / ω
    Y = e*B/mω  # ωₘ²/ω  gyrofrequency / ω
    Z = ν(z)/ω  # collision frequency / ω
    U = 1 - im*Z

    # Temporary variables
    U² = U^2
    Y² = Y^2

    YU = Y*U
    nYU = n*YU
    mYU = m*YU
    lYU = l*YU

    nY² = n*Y²
    lnY² = l*nY²
    mnY² = m*nY²
    lmY² = l*m*Y²

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

    earthcurvature = 2/Rₑ*(H - z)
    M -= earthcurvature*I  # This correction occurs once after summing species

    return M
end

"""
    tmatrix(ea, M)

Return the intermediate matrix components of `T` whose elements are used in the calculation
of matrix `S` for the differential equations for the ionospheric reflection coefficient `R`.

Following Budden's formalism [^Budden1955a] the differential equations for calculating the
reflection coefficient for a radio wave obliquely incident on a horizontally stratified
ionosphere are derived from Maxwell's equations in conjunction with the constitutive relations
for the ionosphere, represented by the matrix `M`. After eliminating the vertically directed
components ``Ez`` and ``Hz``, the equations can be written ``∂e/∂s = -i T e`` where ``s`` is
a proxy for height `z`, and ``e = (Ex -Ey Hx Hy)``, as detailed by [^Clemmow1954]. This
function calculates components of the matrix `T`.

See also: [`smatrix`](@ref), [`susceptibility`](@ref)

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
[^Clemmow1954]: P. C. Clemmow and J. Heading, “Coupled forms of the differential equations governing radio propagation in the ionosphere,” Mathematical Proceedings of the Cambridge Philosophical Society, vol. 50, no. 2, pp. 319–333, Apr. 1954.
"""
function tmatrix(ea::EigenAngle{<:Number}, M)
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
    smatrix(ea, T)

Compute submatrix elements of `S` for solving the differential equation of reflection matrix
`R` wrt `z`.

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
function smatrix(ea::EigenAngle{<:Number}, T::TMatrix{<:Number})
    C = ea.cosθ
    Cinv = inv(C)

    # Temporary matrix elements `T`
    @unpack T11, T12, T14, T31, T32, T34, T41, T42, T44 = T

    # based on Sheddy1968
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

    return S11, S21, S12, S22

    #== The above is equivalent to (but faster than):
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
    dsmatrixdθ(ea, M, T)

Return the 4 submatrix elements of the derivative of the `S` matrix wrt θ.

See also: [`smatrix`](@ref), [`susceptibility`](@ref), [`tmatrix`](@ref)
"""
function dsmatrixdθ(ea::EigenAngle{<:Number}, M, T::TMatrix{<:Number})
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = inv(C)
    C²inv = inv(C²)

    dS = C
    dC = -S
    dC² = -2*S*C
    dCinv = S*C²inv

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

    return dS11, dS21, dS12, dS22
end

"""
    dRdz(R, params, z)

Return the differential of the reflection matrix `R` wrt height `z`.

Following the Budden [^Budden1955a] formalism for the reflection of an (obliquely) incident
plane wave from a horizontally stratified ionosphere, the differential of the reflection
matrix `R` with height `z` can be described by ``2i R′ = S⁽²¹⁾ + S⁽²²⁾R - RS⁽¹¹⁾ - RS⁽¹²⁾R``.
Integrating ``R′`` wrt height `z`, gives the reflection matrix ``R`` for the ionosphere.

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
"""
function dRdz(R, params, z)
    ω, k = params.ω, params.k
    ea = params.ea
    species, bfield = params.species, params.bfield

    M = susceptibility(z, ω, bfield, species)
    T = tmatrix(ea, M)
    S = smatrix(ea, T)

    # the factor k/(2i) isn't explicitly in Budden1955a b/c of his change of variable ``s = kz``
    return -im/2*k*(S[2] + S[4]*R - R*S[1] - R*S[3]*R)
end

function dRdθdz(RdRdθ, params, z)
    ω, k = params.ω, params.k
    ea = params.ea
    species, bfield = params.species, params.bfield

    M = susceptibility(z, ω, bfield, species)
    T = tmatrix(ea, M)
    S = smatrix(ea, T)
    dS = dsmatrixdθ(ea, M, T)

    R = RdRdθ[SVector(1,2),:]
    dRdθ = RdRdθ[SVector(3,4),:]

    dz = -im/2*k*(S[2] + S[4]*R - R*S[1] - R*S[3]*R)
    dθdz = -im/2*k*(dS[2] + dS[4]*R + S[4]*dRdθ -
        (dRdθ*S[1] + R*dS[1]) -
        (dRdθ*S[3]*R + R*dS[3]*R + R*S[3]*dRdθ))

    return vcat(dz, dθdz)
end

function integratedreflection(ea::EigenAngle{<:Number}, integrationparams::IntegrationParameters)
    @unpack bottomheight, topheight, bfield, species, freq = integrationparams

    Mtop = susceptibility(topheight, freq, bfield, species)

    ω, k = freq.ω, freq.k
    params = (ω=ω, k=k, ea=ea, species=species, bfield=bfield)

    Rtop = sharplyboundedreflection(ea, Mtop)

    prob = ODEProblem{false}(dRdz, Rtop, (topheight, bottomheight), params)
    sol = solve(prob, Tsit5(), reltol=1e-6, save_everystep=false, save_start=false)

    return sol[end]
end

function integratedreflection(ea::EigenAngle{<:Number}, integrationparams::IntegrationParameters, ::Derivative)
    @unpack bottomheight, topheight, bfield, species, freq = integrationparams

    Mtop = susceptibility(topheight, freq, bfield, species)

    ω, k = freq.ω, freq.k
    params = (ω=ω, k=k, ea=ea, species=species, bfield=bfield)

    RdRdθtop = sharplyboundedreflection(ea, Mtop, Derivative())

    prob = ODEProblem{false}(dRdθdz, RdRdθtop, (topheight, bottomheight), params)
    sol = solve(prob, Tsit5(), reltol=1e-6, save_everystep=false, save_start=false)

    return sol[end]
end


##########
# Ground reflection coefficient matrix
##########

"""
    fresnelreflection(ea, ground, frequency)

Return the Fresnel reflection coefficient matrix for the ground free-space interface at the
ground (``z = 0``). Follows the formulation in [^Morfitt1976] pages 25-26.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function fresnelreflection(ea::EigenAngle{<:Number}, ground::Ground, frequency::Frequency)
    C, S² = ea.cosθ, ea.sin²θ
    ω = frequency.ω

    Ng² = ground.ϵᵣ - im*ground.σ/(ω*ϵ₀)

    CNg² = C*Ng²
    sqrtNg²S² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²S²)/(CNg² + sqrtNg²S²)
    Rg22 = (C - sqrtNg²S²)/(C + sqrtNg²S²)

    # TODO: Custom type
    return SMatrix{2,2}(Rg11, 0, 0, Rg22)
end

function fresnelreflection(ea::EigenAngle{<:Number}, ground::Ground, frequency::Frequency, ::Derivative)
    C, S, S² = ea.cosθ, ea.sinθ, ea.sin²θ
    ω = frequency.ω

    Ng² = ground.ϵᵣ - im*ground.σ/(ω*ϵ₀)

    CNg² = C*Ng²
    sqrtNg²S² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²S²)/(CNg² + sqrtNg²S²)
    Rg22 = (C - sqrtNg²S²)/(C + sqrtNg²S²)

    dRg11 = (2*Ng²*(1 - Ng²)*S)/(sqrtNg²S²*(CNg² + sqrtNg²S²)^2)
    dRg22 = (2*(C - sqrtNg²S²)*S)/(sqrtNg²S²*(tmsqrtNg²S²p2 + C))

    return SMatrix{2,2}(Rg11, 0, 0, Rg22), SMatrix{2,2}(dRg11, 0, 0, dRg22)
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
    return det(A)*tr(inv(A)*dA)
end

function solvemodalequation(ea::EigenAngle{<:Number}, integrationparams::IntegrationParameters)
    @unpack frequency, ground = integrationparams

    R = integratedreflection(ea, integrationparams)
    Rg = fresnelreflection(ea, ground, frequency)

    f = modalequation(R, Rg)
    return f
end

"""
This returns R and Rg in addition to df because the only time this function is needed, we also
need R and Rg (in excitationfactors).
"""
function solvemodalequationdθ(ea::EigenAngle{<:Number}, integrationparams::IntegrationParameters)
    @unpack frequency, ground = integrationparams

    RdR = integratedreflection(ea, integrationparams, Derivative())
    R = RdR[SVector(1,2),:]
    dR = RdR[SVector(3,4),:]

    Rg, dRg = fresnelreflection(ea, ground, frequency, Derivative())

    df = modalequationdθ(R, dR, Rg, dRg)
    return df, R, Rg
end

function findmodes(origcoords, integrationparams::IntegrationParameters, tolerance=1e-4)
    zroots, zpoles = grpf(z->solvemodalequation(z, integrationparams),
                          origcoords, tolerance, 15000)

    # TODO: do this with SVector?
    return [EigenAngle(r) for r in zroots]
end


################
# Mode Sum Electric Fields
################

"""
MS 76 (pg 38)
Pappert et al 1970 A FORTRAN program...
Pappert and Shockey 1971 WKB Mode Summing Program...

Morfitt 1980 Simplified Eqs 16-23 <=

Assumes `d` = 0.
"""
function heightgainconstants(ea::EigenAngle{<:Number}, ground::Ground, freq::Frequency)
    C², S² = ea.cos²θ, ea.sin²θ
    k, ω = freq.k, freq.ω

    α = 2/Rₑ
    αk23 = (α/k)^(2/3)

    n₀² = 1 - α*H
    Ng² = ground.ϵᵣ - im*ground.σ/(ω*ϵ₀)

    q₀ = (C² - α*H)/αk23

    h₁q₀, h₂q₀, h₁′q₀, h₂′q₀ = modifiedhankel(q₀)

    H₁q₀ = h₁′q₀ + αk23*h₁q₀/2
    H₂q₀ = h₂′q₀ + αk23*h₂q₀/2

    n₀²Ng² = n₀²/Ng²
    cbrtkα_sqrtNg²S² = cbrt(k/α)*sqrt(Ng²-S²)

    F₁ = -(H₂q₀ - im*n₀²Ng²*cbrtkα_sqrtNg²S²*h₂q₀)
    F₂ = H₁q₀ - im*n₀²Ng²*cbrtkα_sqrtNg²S²*h₁q₀
    F₃ = -(h₂′q₀ - im*cbrtkα_sqrtNg²S²*h₂q₀)
    F₄ = h₁′q₀ - im*cbrtkα_sqrtNg²S²*h₁q₀  # NOTE: Typo in P&S 71 eq 9

    return HeightGainConstants(h₁q₀, h₂q₀, h₁′q₀, h₂′q₀, F₁, F₂, F₃, F₄)
end

"""
See Morfitt 1980 Simplified ... Eqs 26-32

Morfitt specifies that `S` is the sine of θ at reference height `h`. The `R`s and `Rg`s
are elements of the reflection matrix looking into iono and towards ground from the level `d`.
Therefore `T`s are also at `d`. BUT (on page 20) "the mode conversion program, as they are
now programmed, require that the height variable `d` be at the ground so that in equations
(26) through (32), `z = d = 0`".

TODO: Return all field components, but have separate path for vertical only.

sw_wvgd.for
"""
function excitationfactors(
    ea::EigenAngle{<:Number},
    R,
    Rg,
    dFdθ,
    hgc::HeightGainConstants,
    fc::FieldComponent
)
    θ, S, S² = ea.θ, ea.sinθ, ea.sin²θ
    sqrtS = sqrt(S)

    @unpack h₁q₀, h₂q₀, h₁′q₀, h₂′q₀, F₁, F₂, F₃, F₄ = hgc

    # Referenced to ground, but `θ` at `h`
    F₁h₁q₀ = F₁*h₁q₀
    F₂h₂q₀ = F₂*h₂q₀
    F₃h₁q₀ = F₃*h₁q₀
    F₄h₂q₀ = F₄*h₂q₀
    D11 = (F₁h₁q₀ + F₂h₂q₀)^2
    D12 = (F₁h₁q₀ + F₂h₂q₀)*(F₃h₁q₀ + F₄h₂q₀)
    D22 = (F₃h₁q₀ + F₄h₂q₀)^2

    # Referenced to ground, but `θ` at `h`
    T₁ = (sqrtS*(1 + Rg[1,1])^2*(1 - R[2,2]*Rg[2,2])) / (dFdθ*Rg[1,1]*D11)
    T₂ = (sqrtS*(1 + Rg[2,2])^2*(1 - R[1,1]*Rg[1,1])) / (dFdθ*Rg[2,2]*D22)
    T₃ = (sqrtS*(1 + Rg[1,1])*(1 + Rg[2,2])*R[2,1]) / (dFdθ*D12)
    T₄ = R[1,2]/R[2,1]

    # This is "undoing" the division by D11, D22, and D12 that are in the definition of the
    # `T`s, but otherwise `f` (below) would be very expensive to calculate directly, so it's
    # worth defining the `T`s as they are.
    τ₁ = D11*T₁
    τ₂ = D22*T₂
    τ₃ = D12*T₃

    # `f` = `Ey/Hy` at the ground (Morfitt 1980, eq. 33)
    # LWPC uses the following criteria for choosing `f`
    if abs2(1 - R[1,1]*Rg[1,1]) > abs2(1 - R[2,2]*Rg[2,2])
        f = T₂/(T₃*T₄)
    else
        f = T₃/T₁
    end

    if fc == FC_Ez
        λ₁ = τ₁*S²
        λ₂ = -τ₁*S
        λ₃ = -τ₃*T₄*S/f
    elseif fc == FC_Ex
        λ₁ = τ₁*S
        λ₂ = -τ₁
        λ₃ = -τ₃*T₄/f
    elseif fc == FC_Ey
        λ₁ = -τ₃*S/f
        λ₂ = τ₃/f
        λ₃ = τ₂/f^2
    end

    return λ₁, λ₂, λ₃, f
end

function heightgains(z, ea::EigenAngle{<:Number}, frequency::Frequency, f, hgc::HeightGainConstants)
    C², S² = ea.cos²θ, ea.sin²θ
    k = frequency.k

    α = 2/Rₑ
    αk23 = (α/k)^(2/3)

    @unpack h₁q₀, h₂q₀, h₁′q₀, h₂′q₀, F₁, F₂, F₃, F₄ = hgc

    q = (C² - α*(H-z))/αk23
    h₁q, h₂q, h₁′q, h₂′q = modifiedhankel(q)

    expzRₑ = exp(z/Rₑ)

    # height gain for vertical electric field Ez, f₁
    f₁ = expzRₑ*(F₁*h₁q + F₂*h₂q)/(F₁*h₁q₀ + F₂*h₂q₀)

    # horizontal electric field Ex, f₂
    f₂ = -im*expzRₑ/(Rₑ*k*(F₁*h₁q₀ + F₂*h₂q₀))

    # Ey, normal to plane of propagation, f₃
    f₃ = (F₃*h₁q + F₄*h₂q)*f/(F₃*h₁q₀ + F₄*h₂q₀)

    return f₁, f₂, f₃
end

function heightgains(z, ea::EigenAngle{<:Number}, frequency::Frequency, f, hgc::HeightGainConstants, fc::FieldComponent)
    C², S² = ea.cos²θ, ea.sin²θ
    k = frequency.k

    α = 2/Rₑ
    αk23 = (α/k)^(2/3)

    @unpack h₁q₀, h₂q₀, h₁′q₀, h₂′q₀, F₁, F₂, F₃, F₄ = hgc

    q = (C² - α*(H-z))/αk23
    h₁q, h₂q, h₁′q, h₂′q = modifiedhankel(q)

    expzRₑ = exp(z/Rₑ)

    if fc == FC_Ez
        # height gain for vertical electric field Ez, f₁
        f = expzRₑ*(F₁*h₁q + F₂*h₂q)/(F₁*h₁q₀ + F₂*h₂q₀)
    elseif fc == FC_Ex
        # horizontal electric field Ex, f₂
        f = -im*expzRₑ/(Rₑ*k*(F₁*h₁q₀ + F₂*h₂q₀))
    elseif fc == FC_Ey
        # Ey, normal to plane of propagation, f₃
        f = (F₃*h₁q + F₄*h₂q)*f/(F₃*h₁q₀ + F₄*h₂q₀)
    end

    return f
end

"""
Morfitt 1980

NOTE: LWPC and Morfitt 1980 have different conventions for (1,2,3). Morfitt has (1,2,3) → (z, x, y)
LWPC has (1,2,3) → (z, y, x) at least for heightgains...

E field strength amplitude in dB above μV/m for 1 kW radiator. Phase relative to free space

see lw_sum_modes.for

amp in dB relative to 1 μV/m for 1000 kW radiator? according to Morfitt 1980, but can't confirm
phase is "relative" phase
"""
function Efield(
    x,
    modes::AbstractVector{EigenAngle{<:Number}},
    integrationparams::IntegrationParameters,
    tx::Exciter,
    rx::AbstractSampler
)

    # Excitation factors require `z = d = 0`
    if integrationparams.bottomheight != 0
        integrationparams = IntegrationParameters(zero(integrationparams.topheight),
                                                  integrationparams.topheight,
                                                  integrationparams.bfield,
                                                  integrationparams.freq,
                                                  integrationparams.ground,
                                                  integrationparams.species)
    end

    @unpack frequency, ground = integrationparams

    talt = altitude(tx)
    ralt = altitude(rx)
    rxfc = fieldcomponent(rx)

    # Transmit dipole antenna orientation with respect to propagation direction
    # See Morfitt 1980 pg 22
    Sγ, Cγ = sincos(elevation(tx))
    Sϕ, Cϕ = sincos(azimuth(tx))

    modesum = zero(ComplexF64)
    for ea in modes
        S = ea.sinθ  # Morfitt 1980 pg 19 specifies S is the sine of θ at `H`

        hgc = heightgainconstants(ea, ground, tx)
        dFdθ, R, Rg = solvemodalequationdθ(ea, integrationparams)

        λ₁, λ₂, λ₃, f = excitationfactors(ea, R, Rg, dFdθ, hgc, rxfc)
        f₁, f₂, f₃ = heightgains(talt, ea, frequency, f, hgc)  # at transmitter

        # All 3 directions needed for xmtr, only "wanted" term needed for rcvr
        xmtrterm = λ₁*f₁*Cγ + λ₂*f₂*Sγ*Cϕ + λ₃*f₃*Sγ*Sϕ
        rcvrterm = heightgains(ralt, ea, frequency, f, hgc, rxfc)  # at receiver

        # NOTE: Morfitt 1980 eq 39 is different from Pappert et al 1983 which is also
        # different from Pappert and Ferguson 1986
        modesum += xmtrterm*rcvrterm*cis(-k*(S - 1)*x)
    end

    # Q = 3.248e-5*sqrt(1000)*k/sqrt(f)  # Morfitt 1980 eq 41, adjusted for MKS units. power is also probably included in this
    # Express field components in volts/meter/kW of radiated power
    Q = 6.808e-4*sqrt(f/1000)  # Pappert and Ferguson 1986, pg 553, adjusted for MKS

    # Total electric field
    E = Q/sqrt(sin(x/Rₑ))*modesum
    phase = angle(E)
    amp = 10*log10(norm(E))

    return E, phase, amp
end

"""
Special function for vertical Ez only, zt=zr=0, and γ=φ=0. (Like GRNDMC)
essentially, zt=zr=0, f₁=f₂=f₃=1, γ=ϕ=0

TODO: Efield function wrapper should automatically call this one if criteria are met

Q is also in Pappert Snyder 1972
"""
function Efield(x, modes::AbstractVector{EigenAngle{<:Number}}, integrationparams::IntegrationParameters)
    @unpack tx, ground = integrationparams
    k, f = tx.k, tx.freq

    # Q = 3.248e-5*sqrt(1000)*k/sqrt(f)  # Morfitt 1980 eq 41, adjusted for MKS units. power is also probably included in this
    # QM = 9.023e-8*complex(f)^(3/2)  #  Pappert et al 1983
    # Q = 0.03248*k/sqrt(f)*sqrt(1000)  # TODO: find an explicit reference for this
    Q = 6.808e-4*sqrt(f/1000)  # Pappert and Ferguson 1986, pg 553

    modesum = zero(ComplexF64)
    for eaₕ in modes
        ea₀ = referenceground(eaₕ)
        # ea₀ = eaₕ # BUG: why does this work better than including the conversion???
        S = ea₀.sinθ

        hgc = heightgainconstants(eaₕ, ground, tx)

        # Excitation factors require `z = d = 0`
        if integrationparams.bottomheight != 0
            integrationparams = IntegrationParameters(zero(integrationparams.topheight),
                                                      integrationparams.topheight,
                                                      integrationparams.bfield,
                                                      integrationparams.tx,
                                                      integrationparams.ground,
                                                      integrationparams.species)
        end
        dFdθ, R, Rg = solvemodalequationdθ(eaₕ.θ, integrationparams)

        λ_Ez, λ_Ex, λ_Ey = excitationfactors(eaₕ, R, Rg, dFdθ, hgc, Ez())

        xmtrterm = λ_Ez

        modesum += xmtrterm*cis(-k*(S - 1)*x)
    end

    E = Q/sqrt(abs(sin(x/Rₑ)))*modesum
    phase = angle(E)
    amp = 10*log10(norm(E))

    return E, phase, amp
end

function Efield(
    X::AbstractArray,
    modes::AbstractVector{EigenAngle{<:Number}},
    integrationparams::IntegrationParameters
)
    Xlength = length(X)
    E = Array{ComplexF64}(undef, Xlength)
    phase = Array{Float64}(undef, Xlength)
    amp = Array{Float64}(undef, Xlength)
    for i in eachindex(X)
        e, p, a = Efield(X[i], modes, integrationparams)
        E[i] = e
        phase[i] = p  # TODO: Wrap phase and find explicit reference (e.g. wrt antenna current phase)
        amp[i] = a
    end
    return E, phase, amp
end

function Efield(
    X::AbstractArray,
    modes::AbstractVector{EigenAngle{<:Number}},
    integrationparams::IntegrationParameters,
    rx::AbstractReceiver,
    fc::FieldComponent
)
    Xlength = length(X)
    E = Array{ComplexF64}(undef, Xlength)
    phase = Array{Float64}(undef, Xlength)
    amp = Array{Float64}(undef, Xlength)
    for i in eachindex(X)
        e, p, a = Efield(X[i], modes, integrationparams, rx, fc)
        E[i] = e
        phase[i] = p  # TODO: Wrap phase and find explicit reference (e.g. wrt antenna current phase)
        amp[i] = a
    end
    return E, phase, amp
end

"""
Convert `ea` at ground to `ea` at reference height `H`.

See Morfitt 1980 pg 26
"""
function referenceheight(ea::EigenAngle{<:Number})
    θₕ = asin((1 - H/Rₑ)*ea.sinθ)
    return EigenAngle(θₕ)
end

function referenceground(ea::EigenAngle{<:Number})
    θ₀ = asin(ea.sinθ/(1 - H/Rₑ))
    return EigenAngle(θ₀)
end

"""
Sheddy 1968 specifies this is the phase velocity at the ground (ea at ground)

MS76 pg 82

Sheddy includes K where K = 1 + α H/2 presumably to correct for earth curvature

See Pappert et al 1967 Numerical Investigation...

in m/s
"""
function phasevelocity(ea::EigenAngle{<:Number})
    K = 1 + H/Rₑ
    return c₀*K/real(ea.sinθ)
end

"""
ms 76 pg 82 and sheddy 1968

Sheddy includes K where K = 1 + α H/2 presumably to correct for earth curvature

The strange factor in front converts 1 neper to dB

See Pappert et al 1967 Numerical Investigation...


in dB/Mm
"""
function attenuationrate(ea::EigenAngle{<:Number}, tx::AbstractSource)
    K = 1 + H/Rₑ
    return -20*log10(exp(1))*K*tx.k*imag(ea.sinθ)
end
