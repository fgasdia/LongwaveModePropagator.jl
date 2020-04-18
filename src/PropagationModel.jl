const TOPHEIGHT = 95e3  # temporary - should be part of an actual IntegrationParameters
const BOTTOMHEIGHT = zero(TOPHEIGHT)  # should this be an actual const? Nothing works if it's not 0...

#==
A struct like this could be used in place of the `const`s below.
It would allow us to parameterize the Complex type, but realistically it will
never be anything other than ComplexF64.

# PolynomialRoots package requires complex floats of arbitrary precision
struct BookerQuartic{T<:Complex{<:AbstractFloat}}
    roots::MVector{4,T}
    coeffs::MVector{5,T}
end
==#

# Passing MArrays between functions causes allocations. They are avoided by
# mutating this const in place. `roots!` requires Complex values.
const BOOKER_QUARTIC_ROOTS = MVector{4}(zeros(ComplexF64, 4))
const BOOKER_QUARTIC_COEFFS = MVector{5,ComplexF64}(undef)

struct Derivative_dθ end

@with_kw struct ModeParameters{F, G}
    bfield::BField
    frequency::Frequency
    ground::Ground
    species::Constituent{F,G}
end

@with_kw struct FresnelReflectionVariables{T<:Number} @deftype T
    Ng²
    CNg²
    sqrtNg²mS²
    Rg::SMatrix{2,2,T,4}
end

# TODO: G as its own type contained within this?
@with_kw struct SharpBoundaryVariables{T<:Number} @deftype T
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

@with_kw struct ExcitationFactor{T} @deftype T
    F₁
    F₂
    F₃
    F₄
    h₁0
    h₂0
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
function bookerquartic!(ea::EigenAngle, M::AbstractArray{ComplexF64})
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # Precompute
    M11p1 = 1 + M[1,1]
    M33p1 = 1 + M[3,3]

    C²pM22 = C² + M[2,2]
    C²pM33 = C² + M[3,3]

    M13pM31 = M[1,3] + M[3,1]

    M12M23 = M[1,2]*M[2,3]
    M13M31 = M[1,3]*M[3,1]
    M21M32 = M[2,1]*M[3,2]
    M23M32 = M[2,3]*M[3,2]

    # Booker quartic coefficients
    B4 = M33p1
    B3 = S*M13pM31
    B2 = -C²pM33*M11p1 + M13M31 -
         M33p1*C²pM22 + M23M32
    B1 = S*(M12M23 + M21M32 -
         C²pM22*M13pM31)
    B0 = M11p1*C²pM22*C²pM33 +
         M12M23*M[3,1] + M[1,3]*M21M32 -
         M13M31*C²pM22 -
         M11p1*M23M32 -
         M[1,2]*M[2,1]*C²pM33

    BOOKER_QUARTIC_COEFFS.data = (B0, B1, B2, B3, B4)

    # NOTE: `roots!` dominates this functions runtime
    # It may be worthwhile to replace this with my own or improve PolynomialRoots.
    roots!(BOOKER_QUARTIC_ROOTS, BOOKER_QUARTIC_COEFFS, NaN, 4, false)

    return nothing
end

"""
    bookerquartic(T::TMatrix)

Solves the depressed quartic in terms of `T`. This technique is used in LWPC's
wavefield subroutines, e.g. "wf_init.for".

This function is ~2× as fast as [`bookerquartic(ea, M)`](@ref).
"""
function bookerquartic!(T::TMatrix{ComplexF64})
    # This is the depressed form of the quartic
    b3 = -(T[1,1] + T[4,4])
    b2 = T[1,1]*T[4,4] - T[1,4]*T[4,1] - T[3,2]
    b1 = -(-T[3,2]*(T[1,1] + T[4,4]) + T[1,2]*T[3,1] + T[3,4]*T[4,2])
    b0 = -T[1,1]*(T[3,2]*T[4,4] - T[3,4]*T[4,2]) +
        T[1,2]*(T[3,1]*T[4,4] - T[3,4]*T[4,1]) -
        T[1,4]*(T[3,1]*T[4,2] - T[3,2]*T[4,1])

    BOOKER_QUARTIC_COEFFS.data = (b0, b1, b2, b3, one(ComplexF64))

    roots!(BOOKER_QUARTIC_ROOTS, BOOKER_QUARTIC_COEFFS, NaN, 4, false)

    return nothing
end


"""
    upgoing(v)

Calculate the absolute angle of `v` in radians from 315°×π/180 on the
complex plane. Smaller values indicate upgoing waves.
"""
function upgoing(v::Complex)
    # -π/4 is 315°
    a = -π/4 - angle(v)

    # angle() is +/-π
    # `a` can never be > 3π/4, but there is a region on the plane where `a` can be < -π
    a < -π && (a += 2π)
    return abs(a)
end

function upgoing(v::Real)
    return oftype(v, π/4)
end

"""

          Im
      3 . |
          |
          |
  . 4     |
 ------------------Re
          |      . 2
          |
          |
          | . 1

General locations of the four roots of the Booker quartic `q`:
    1) upgoing evanescent wave
    2) upgoing travelling wave
    3) downgoing evanescent wave
    4) downgoing travelling wave

From Pitteway 1967 fig 5. Sheddy 1968 also mentions that the two upgoing waves
correspond to the two solutions that lie close to the positive real axis and the
negative imaginary axis, respectively.

NOTE: It looks like we could just `sort!(q, by=real)`, but the quadrants of each
root are not fixed and the positions in the Argand diagram are just approximate.

See also [`sortquarticroots!`](@ref)
"""
function sortquarticroots(v::AbstractArray{T}) where {T <: Complex}
    # Array where we order roots
    swap = copy(v)

    sortquarticroots!(swap)

    if v isa MVector
        return SVector(swap)
    else
        return swap
    end
end

function sortquarticroots!(v::AbstractArray{Complex{T}}) where {T}

    # Calculate and sort by distance from 315°
    # The two closest are upgoing and the two others are downgoing
    # This is faster than `sort!(v, by=upgoing)`
    dist = MVector{length(v),T}(undef)
    for i in eachindex(v)
        dist[i] = upgoing(v[i])
    end

    i = length(v)
    while i > 1
        if dist[i] < dist[i-1]
            dist[i], dist[i-1] = dist[i-1], dist[i]
            v[i], v[i-1] = v[i-1], v[i]
            i = length(v) # Need to restart at the end
        else
            i -= 1
        end
    end

    # Now arrange the roots so that they correspond to:
    # 1) upgoing evanescent, 2) upgoing travelling
    # 3) downgoing evanescent, and 4) downgoing travelling
    if imag(v[1]) > imag(v[2])
        v[2], v[1] = v[1], v[2]
    end
    if imag(v[3]) < imag(v[4])
        v[4], v[3] = v[3], v[4]
    end

    return v
end


function _sharpboundaryreflection(ea::EigenAngle, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # XXX: `bookerquartic` (really `roots!`) dominates this functions runtime
    bookerquartic!(ea, M)

    # We choose the 2 roots corresponding to upward travelling waves as being those that lie
    # close to the positive real and negative imaginary axis (315° on the complex plane)
    sortquarticroots!(BOOKER_QUARTIC_ROOTS)

    #==
    Γ = [0 -q 0;
         q 0 -S;
         0 S 0]
    G = Γ² + I + M
    G = [1-q[i]^2+M[1,1] M[1,2] S*q[i]+M[1,3];
         M[2,1] 1-q[i]^2-S²+M[2,2] M[2,3];
         S*q[i]+M[3,1] M[3,2] C²+M[3,3]]
    ==#

    #BUG: Need to check that Sheddy actually wants roots 2, 1 in that order rather
    # than the order used by Pitteway.

    # Precompute
    q = BOOKER_QUARTIC_ROOTS
    q₁S = q[1]*S
    q₂S = q[2]*S

    M11p1 = 1 + M[1,1]

    # Constant entries of dispersion matrix `G`
    G12 = M[1,2]
    G32 = M[3,2]
    G33 = C² + M[3,3]

    G12G33 = G12*G33

    # Values for two solutions of the Booker quartic corresponding to upgoing waves
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

    # Computation of entries of reflection matrix `R`
    Δ = (T₁*C + P₁)*(C + q[2]) - (T₂*C + P₂)*(C + q[1])
    Δ⁻¹ = 1/Δ

    return SharpBoundaryVariables(G12, G32, G33, G11₁, G13₁, G31₁, Δ₁, Δ₁⁻¹, P₁, T₁,
                                  G11₂, G13₂, G31₂, Δ₂, Δ₂⁻¹, P₂, T₂, Δ, Δ⁻¹)
end

"""
    sharplyboundedreflection(ea, M)

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

    @unpack P₁, T₁, P₂, T₂, Δ⁻¹ = _sharpboundaryreflection(ea, M)
    q = sortquarticroots(BOOKER_QUARTIC_ROOTS)

    T₁C = T₁*C
    T₂C = T₂*C

    R11 = ((T₁C - P₁)*(C + q[2]) - (T₂C - P₂)*(C + q[1]))*Δ⁻¹  # ∥R∥
    R22 = ((T₁C + P₁)*(C - q[2]) - (T₂C + P₂)*(C - q[1]))*Δ⁻¹  # ⟂R⟂
    R12 = -C2*(T₁*P₂ - T₂*P₁)*Δ⁻¹  # ⟂R∥
    R21 = -C2*(q[1] - q[2])*Δ⁻¹  # ∥R⟂

    return SMatrix{2,2}(R11, R21, R12, R22)
end

# TODO: Autodiff with Zygote?
function sharpboundaryreflection(ea::EigenAngle, M, ::Derivative_dθ)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    @unpack G12, G32, G33, G11₁, G13₁, G31₁, Δ₁, Δ₁⁻¹, P₁, T₁,
        G11₂, G13₂, G31₂, Δ₂, Δ₂⁻¹, P₂, T₂, Δ, Δ⁻¹ = _sharpboundaryreflection(ea, M)
    q, B = BOOKER_QUARTIC_ROOTS, BOOKER_QUARTIC_COEFFS

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

    return @SMatrix [R11 R12;
                     R21 R22;
                     dR11 dR12;
                     dR21 dR22]
end

##########
# Reflection coefficient matrix for a vertically stratified ionosphere
##########

"""
    susceptibility(z, frequency, bfield, species)

Computation of susceptibility matrix `M` as defined by [^Budden1955a] method 2.

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
[^Budden1988]: K. G. Budden
"""
function susceptibility(z, frequency::Frequency, bfield::BField, species::Constituent)
    B, lx, ly, lz = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    lx², ly², lz² = lx^2, ly^2, lz^2
    ω = frequency.ω

    # Constitutive relations (see Budden1955a, pg. 517 or Budden1988 pg. 39)
    e, m, N, ν = species.charge, species.mass, species.numberdensity, species.collisionfrequency
    mω = m*ω
    X = N(z)*e^2/(ϵ₀*mω*ω)  # ωₚ²/ω² plasma frequency / ω
    Y = e*B/mω  # ωₘ²/ω  gyrofrequency / ω  # BUG: What if e is negative?
    Z = ν(z)/ω  # collision frequency / ω
    U = 1 - im*Z

    # TODO: Support multiple species, e.g.
    # U²D = zero()
    # Y²D = zero()
    # YUD = zero()
    # for i = 1:length(species)
    #     U²D += U^2*D
    #     YUD += Y*U*D
    #     Y²D += Y^2*D
    # end

    # Precompute variables
    # These are formulated such that application to multiple species is trivial. Direction
    # cosines are applied last
    U² = U^2
    Y² = Y^2

    D = -X/(U*(U² - Y²))

    U²D = U²*D
    Y²D = Y²*D
    YUD = Y*U*D

    # Leverage partial symmetry of M to reduce computations
    izYUD = im*lz*YUD
    xyY²D = lx*ly*Y²D
    iyYUD = im*ly*YUD
    xzY²D = lx*lz*Y²D
    ixYUD = im*lx*YUD
    yzY²D = ly*lz*Y²D

    earthcurvature = 2/Rₑ*(H - z)

    # Elements of `M`
    M11 = U²D - lx²*Y²D - earthcurvature
    M21 = izYUD - xyY²D
    M31 = -iyYUD - xzY²D
    M12 = -izYUD - xyY²D
    M22 = U²D - ly²*Y²D - earthcurvature
    M32 = ixYUD - yzY²D
    M13 = iyYUD - xzY²D
    M23 = -ixYUD - yzY²D
    M33 = U²D - lz²*Y²D - earthcurvature

    # Remember, column major
    M = SMatrix{3,3}(M11, M21, M31,
                     M12, M22, M32,
                     M13, M23, M33)

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

See also: [`wmatrix`](@ref), [`susceptibility`](@ref)

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
[^Clemmow1954]: P. C. Clemmow and J. Heading, “Coupled forms of the differential equations governing radio propagation in the ionosphere,” Mathematical Proceedings of the Cambridge Philosophical Society, vol. 50, no. 2, pp. 319–333, Apr. 1954.
"""
function tmatrix(ea::EigenAngle, M)
    S, C² = ea.sinθ, ea.cos²θ

    # Denominator of most of the entries of `T`
    # This expression dominates the function runtime
    den = inv(1 + M[3,3])

    M31den = M[3,1]*den
    M32den = M[3,2]*den

    T11 = -S*M31den
    T12 = S*M32den
    # T13 = 0
    T14 = (C² + M[3,3])*den
    # T21 = 0
    # T22 = 0
    # T23 = 1
    # T24 = 0
    T31 = M[2,3]*M31den - M[2,1]
    T32 = C² + M[2,2] - M[2,3]*M32den
    # T33 = 0
    T34 = S*M[2,3]*den
    T41 = 1 + M[1,1] - M[1,3]*M31den
    T42 = M[1,3]*M32den - M[1,2]
    # T43 = 0
    T44 = -S*M[1,3]*den

    # Remember, column major
    return TMatrix(T11, T31, T41,
                   T12, T32, T42,
                   T14, T34, T44)
end

"""
    wmatrix(ea, T)

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
function wmatrix(ea::EigenAngle, T::TMatrix)
    C = ea.cosθ
    Cinv = inv(C)

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
    W11 = @SMatrix [a11+a11r -b11;
                    -c11 d11]
    W12 = @SMatrix [a21+a21r -b11;
                    c12 d12]
    W21 = @SMatrix [a21-a21r b21;
                    c11 -d12]
    W22 = @SMatrix [a11-a11r b21;
                    -c12 -d11]

    return W11, W21, W12, W22
end

"""
    dwmatrixdθ(ea, M, T)

Return the 4 submatrix elements of the derivative of the `W` matrix wrt θ.

See also: [`wmatrix`](@ref), [`susceptibility`](@ref), [`tmatrix`](@ref)
"""
function dwmatrixdθ(ea::EigenAngle, M, T::TMatrix)
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = inv(C)
    C²inv = inv(C²)

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
    dW11 = @SMatrix [ds11a+ds11b -ds12;
                     -ds21 ds22]
    dW12 = @SMatrix [-dd11a+dd11b -ds12;
                     dd21 -dd22]
    dW21 = @SMatrix [-dd11a-dd11b dd12;
                     ds21 dd22]
    dW22 = @SMatrix [ds11a-ds11b dd12;
                     -dd21 -ds22]

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
    ea, modeparams = params
    @unpack bfield, frequency, species = modeparams

    k = frequency.k

    # XXX: `susceptibility` takes as much time as `tmatrix` and `wmatrix` combined
    # TODO: Maybe use a function approximator version of M b/c it's not dependent on θ
    M = susceptibility(z, frequency, bfield, species)
    T = tmatrix(ea, M)
    W = wmatrix(ea, T)

    # the factor k/(2i) isn't explicitly in Budden1955a b/c of his change of variable ``s = kz``
    return -im/2*k*(W[2] + W[4]*R - R*W[1] - R*W[3]*R)
end

function dRdθdz(RdRdθ, params, z)
    ea, modeparams = params
    @unpack bfield, frequency, species = modeparams

    k = frequency.k

    M = susceptibility(z, frequency, bfield, species)
    T = tmatrix(ea, M)
    W = wmatrix(ea, T)
    dW = dwmatrixdθ(ea, M, T)

    R = RdRdθ[SVector(1,2),:]
    dRdθ = RdRdθ[SVector(3,4),:]

    dz = -im/2*k*(W[2] + W[4]*R - R*W[1] - R*W[3]*R)
    dθdz = -im/2*k*(dW[2] + dW[4]*R + W[4]*dRdθ -
            (dRdθ*W[1] + R*dW[1]) -
            (dRdθ*W[3]*R + R*dW[3]*R + R*W[3]*dRdθ))

    return vcat(dz, dθdz)
end

function integratedreflection(ea::EigenAngle, modeparams::ModeParameters)
    @unpack bfield, frequency, species = modeparams

    Mtop = susceptibility(TOPHEIGHT, frequency, bfield, species)

    Rtop = sharpboundaryreflection(ea, Mtop)

    prob = ODEProblem{false}(dRdz, Rtop, (TOPHEIGHT, BOTTOMHEIGHT), (ea, modeparams))
    sol = solve(prob, BS5(), abstol=1e-6, reltol=1e-6, save_everystep=false, save_start=false)

    return sol[end]
end

# This is kept as a completely separate function because the size of the matrix being
# integrated is different and therefore the size of sol[end] is different too
# The derivative terms are intertwined with the non-derivative terms so we can't do only
# the derivative terms
function integratedreflection(ea::EigenAngle, modeparams::ModeParameters, ::Derivative_dθ)
    @unpack bfield, frequency, species = modeparams

    Mtop = susceptibility(TOPHEIGHT, frequency, bfield, species)

    RdRdθtop = sharpboundaryreflection(ea, Mtop, Derivative_dθ())

    prob = ODEProblem{false}(dRdθdz, RdRdθtop, (TOPHEIGHT, BOTTOMHEIGHT), (ea, modeparams))
    sol = solve(prob, BS5(), abstol=1e-6, reltol=1e-6, save_everystep=false, save_start=false)

    return sol[end]
end


##########
# Ground reflection coefficient matrix
##########

function _fresnelreflection(ea, ground, frequency)
    C, S² = ea.cosθ, ea.sin²θ
    ω = frequency.ω

    Ng² = ground.ϵᵣ - im*ground.σ/(ω*ϵ₀)

    CNg² = C*Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²mS²)/(CNg² + sqrtNg²mS²)
    Rg22 = (C - sqrtNg²mS²)/(C + sqrtNg²mS²)

    Rg = SMatrix{2,2}(Rg11, 0, 0, Rg22)

    return FresnelReflectionVariables(Ng², CNg², sqrtNg²mS², Rg)
end

"""
    fresnelreflection(ea, ground, frequency)

Return the Fresnel reflection coefficient matrix for the ground free-space interface at the
ground (``z = 0``). Follows the formulation in [^Morfitt1976] pages 25-26.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency)
    @unpack Rg = _fresnelreflection(ea, ground, frequency)

    return Rg
end

function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency, ::Derivative_dθ)
    C, S, S² = ea.cosθ, ea.sinθ, ea.sin²θ
    ω = frequency.ω

    @unpack Rg, Ng², CNg², sqrtNg²mS² = _fresnelreflection(ea, ground, frequency)

    S2 = 2*S

    dRg11 = (S2*Ng²*(1 - Ng²))/(sqrtNg²mS²*(CNg² + sqrtNg²mS²)^2)
    dRg22 = (S2*(C - sqrtNg²mS²))/(sqrtNg²mS²*(sqrtNg²mS² + C))

    dRg = SMatrix{2,2}(dRg11, 0, 0, dRg22)

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
    return det(A)*tr(inv(A)*dA)
end

function solvemodalequation(ea::EigenAngle, modeparams::ModeParameters)
    @unpack frequency, ground = modeparams

    R = integratedreflection(ea, modeparams)
    Rg = fresnelreflection(ea, ground, frequency)

    f = modalequation(R, Rg)
    return f
end

"""
This returns R and Rg in addition to df because the only time this function is needed, we also
need R and Rg (in excitationfactors).
"""
function solvemodalequationdθ(ea::EigenAngle, modeparams::ModeParameters)
    @unpack frequency, ground = modeparams

    RdR = integratedreflection(ea, modeparams, Derivative_dθ())
    R = RdR[SVector(1,2),:]
    dR = RdR[SVector(3,4),:]

    Rg, dRg = fresnelreflection(ea, ground, frequency, Derivative_dθ())

    df = modalequationdθ(R, dR, Rg, dRg)
    return df, R, Rg
end

"""
`tolerance` is how close solution of modal equation gets to 0. Difference in θ between
`tolerance=1e-6` and `tolerance=1e-8` is less than 1e-5° in both real and imaginary
components.
"""
function findmodes(origcoords, modeparams::ModeParameters, tolerance=1e-6)
    angletype = eltype(origcoords)
    zroots, zpoles = grpf(θ->solvemodalequation(EigenAngle{angletype}(θ), modeparams),
                          origcoords, tolerance, 30000)

    return [EigenAngle{angletype}(r) for r in zroots]
end

################
# Wavefields
################


################
# Excitation and Mode Sum
################

"""
TODO: LF corrections/fixes

Pappert 1981 LF Daytime Earth Ionosphere...
finds that the linear combination of modified Hankel functions of order one third used to
represent the height gain at the ground is in fact incorrect at the ground for modes which
are highly earth detached (pg 7-8). To correct for this, they just throw away earth
effects altogether when the condition
Re(2im*(k/α)*(C²ₕ - α*H)^(3/2)/3)) > 12.4
is met. The value of 12.4 was found by trial and error and requires the degree of evanescence
at the ground to be of the order of several times 10⁻⁶. When this condition is met, the
plane wave reflection coefficients (for ground) become eq (2) and (3).
There are additional equation replacements on page 10.
Otherwise, this is following the math from Pappert & Shockey 71. This paper (81) explicitly
identifies where some of the angles should be referenced.
"""

"""
TODO: ELF corrections/fixes

Pappert Shockey 1971 WKB Mode Summing...
changes height gain functions for ELF, pg. 9

Switches to flat earth at imaginary angles less than -10° (see LWPC or Pappert 1983 ELF-VLF)
"""

"""
Pappert Shockey 1971 pg 9 or Pappert 198e pg 12
"""
function elf_heightgains()

end

"""
    pow23

Calculate `x^(2/3)` relatively efficiently as ``exp(2/3*log(x))``.
"""
pow23(x) = exp(2/3*log(x))

"""
    excitationfactorconstants(ea, R, Rg, frequency, ground)

Return an `ExcitationFactor` struct used in calculating height gains.

Based on Morfitt 1980, Pappert Shockey 1971, and Pappert Shockey 1976 (this last one has H=0)

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Ferguson1981] J. A. Ferguson and D. G. Morfitt, “WKB mode summing program for dipole antennas of arbitrary orientation and elevation for VLF/LF propagation,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-697, Oct. 1981.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function excitationfactorconstants(ea::EigenAngle, R, Rg, frequency::Frequency, ground::Ground)
    S², C² = ea.sin²θ, ea.cos²θ
    k, ω = frequency.k, frequency.ω

    # `ea` is at height `H`. See, e.g. Pappert1971 pg 8

    # Precompute
    α = 2/Rₑ
    αH = α*H
    koα = k/α
    cbrtkoα = cbrt(koα)
    koα23 = cbrtkoα^2  # (k/α)^(2/3)
    αok23 = inv(koα23)  # (α/k)^(2/3)

    q₀ = koα23*(C² - αH)  # (k/α)^(2/3)*(C² - αH)

    # XXX: `modifiedhankel` dominates runtime of this function
    h₁0, h₂0, h₁p0, h₂p0 = modifiedhankel(q₀)

    H₁0 = h₁p0 + αok23*h₁0/2
    H₂0 = h₂p0 + αok23*h₂0/2

    n₀² = 1 - αH  # modified index of refraction (free space) squared
    Ng² = ground.ϵᵣ - im*ground.σ/(ω*ϵ₀)  # ground index of refraction

    # Precompute
    n₀²oNg² = n₀²/Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    cbrtkoαsqrtNg²mS²h₂0 = cbrtkoα*sqrtNg²mS²*h₂0
    cbrtkoαsqrtNg²mS²h₁0 = cbrtkoα*sqrtNg²mS²*h₁0

    F₁ = -H₂0 + im*n₀²oNg²*cbrtkoαsqrtNg²mS²h₂0
    F₂ = H₁0 - im*n₀²oNg²*cbrtkoαsqrtNg²mS²h₁0
    F₃ = -h₂p0 + im*cbrtkoαsqrtNg²mS²h₂0
    F₄ = h₁p0 - im*cbrtkoαsqrtNg²mS²h₁0

    # ey/hy; polarization ratio; Normalizes y component of H to unity at thr ground.
    # Sometimes called `f` in papers
    # f0fr = T₂/(T₃*T₄) = T₃/T₁
    # LWPC uses the following criteria for choosing
    # if abs2(1 - R[1,1]*Rg[1,1]) > abs2(1 - R[2,2]*Rg[2,2])
    #     f0fr = (1 + Rg[2,2])*(1 - R[1,1]*Rg[1,1])/((1 + Rg[1,1])*R[1,2]*Rg[2,2])
    # else
    #     f0fr = (1 + Rg[2,2])*R[2,1]*Rg[1,1]/((1 + Rg[1,1])*(1 - R[2,2]*Rg[2,2]))
    # end

    return ExcitationFactor(F₁, F₂, F₃, F₄, h₁0, h₂0)
end

"""
    heightgains(z, ea, frequency, efconstants)

Calculate heightgains at `z`.

This function assumes the reference height for the reflection coefficients is ``d=0``.

See also: [`excitationfactor`](@ref), [`excitationfactorconstants`](@ref)

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function heightgains(z, ea::EigenAngle, frequency::Frequency, efconstants::ExcitationFactor)
    C² = ea.cos²θ
    k = frequency.k
    @unpack F₁, F₂, F₃, F₄, h₁0, h₂0 = efconstants

    # Precompute
    α = 2/Rₑ
    koα23 = pow23(k/α)  # (k/α)^(2/3)
    expzoRₑ = exp(z/Rₑ)  # assumes `d = 0`

    qz = koα23*(C² - α*(H - z))

    # XXX: `modifiedhankel` dominates this functions runtime
    h₁z, h₂z, h₁pz, h₂pz = modifiedhankel(qz)

    # Precompute
    F₁h₁0 = F₁*h₁0
    F₂h₂0 = F₂*h₂0
    F₁h₁z = F₁*h₁z
    F₂h₂z = F₂*h₂z

    # From Ferguson1981 pg 10 or Pappert1971 pg 6, 8:
    # Height gain for Ez, also called f∥(z)
    fz = expzoRₑ*(F₁h₁z + F₂h₂z)

    # Height gain for Ex, also called g(z)
    # f₂ = 1/(im*k) df₁/dz
    fx = -im*expzoRₑ/(Rₑ*k)*(F₁h₁z + F₂h₂z + Rₑ*(F₁*h₁pz + F₂*h₂pz))

    # Height gain for Ey, also called f⟂(z)
    fy = (F₃*h₁z + F₄*h₂z)

    return fz, fx, fy
end

"""
    excitationfactor(ea, dFdθ, R, Rg, component)

Calculate the excitation factor for electric field `component`.

These excitation factors are used in conjunction with the function [`heightgains`](@ref).

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Ferguson1981] J. A. Ferguson and D. G. Morfitt, “WKB mode summing program for dipole antennas of arbitrary orientation and elevation for VLF/LF propagation,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-697, Oct. 1981.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function excitationfactor(ea::EigenAngle, dFdθ, R, Rg, efconstants::ExcitationFactor, component::FieldComponent)
    S, S², C² = ea.sinθ, ea.sin²θ, ea.cos²θ
    sqrtS = sqrt(S)

    @unpack F₁, F₂, F₃, F₄, h₁0, h₂0 = efconstants

    # TODO: Beter incorporation with other functions?
    # Calculate height gains for `z=d=0`
    fz = (F₁*h₁0 + F₂*h₂0)
    fy = (F₃*h₁0 + F₄*h₂0)

    D₁₁ = fz^2
    D₁₂ = fz*fy
    D₂₂ = fy^2

    # `S` is at `H, specified in e.g.
    # D. G. Morfitt, ``'Simplified' VLF/LF mode conversion...,'' NOSC/TR-514, 1980, pg 19.
    T₁ = sqrtS*(1 + Rg[1,1])^2*(1 - R[2,2]*Rg[2,2])/(dFdθ*Rg[1,1]*D₁₁)
    T₂ = sqrtS*(1 + Rg[2,2])^2*(1 - R[1,1]*Rg[1,1])/(dFdθ*Rg[2,2]*D₂₂)
    T₃ = sqrtS*(1 + Rg[1,1])*(1 + Rg[2,2])*R[2,1]/(dFdθ*D₁₂)
    T₄ = R[1,2]/R[2,1]

    ST₁ = S*T₁
    ST₃ = S*T₃
    if component == FC_Ez
        λv = S²*T₁
        λb = -ST₃*T₄
        λe = -ST₁
    elseif component == FC_Ey
        λv = -ST₃
        λb = T₂
        λe = T₃
    elseif component == FC_Ex
        λv = ST₁
        λb = -T₃*T₄
        λe = -T₁
    end

    return λv, λb, λe
end

"""
    Efield(x, modes, modeparams, tx, rx)

Calculate the complex electric field, amplitude, and phase at a distance `x` from transmitter `tx`.

`modes` is a collection of `EigenAngles` for the earth-ionosphere waveguide with the parameters
`modeparams`. Exciter `tx` specifies the transmitting antenna position, orientation, and
radiated power, and `rx` specifies the field component of interest.

# References


"""
function Efield(x, modes, modeparams::ModeParameters, tx::Exciter, rx::AbstractSampler)
    @unpack frequency, ground = modeparams

    # TODO: special function for vertical component, transmitter, and at ground

    txpower = power(tx)
    zₜ = altitude(tx)
    k = frequency.k
    zᵣ = altitude(rx)
    rxcomponent = fieldcomponent(rx)

    # Transmit dipole antenna orientation with respect to propagation direction
    # See Morfitt 1980 pg 22
    # TODO: Confirm alignment of coord system and magnetic field
    Sγ, Cγ = sincos(π/2 - elevation(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(azimuth(tx))  # ϕ is measured from `x`

    modesum = zero(eltype(eltype(modes)))  # probably ComplexF64
    for ea in modes
        S₀ = ea.sinθ/sqrt(1 - 2/Rₑ*H)  # S referenced to ground, see, e.g. PS71 pg 11

        dFdθ, R, Rg = solvemodalequationdθ(ea, modeparams)
        efconstants = excitationfactorconstants(ea, R, Rg, frequency, ground)

        # fz0, fx0, fy0 = heightgains(0, ea, frequency, efconstants)
        λv, λb, λe = excitationfactor(ea, dFdθ, R, Rg, efconstants, rxcomponent)

        # Transmitter term
        fzₜ, fxₜ, fyₜ = heightgains(zₜ, ea, frequency, efconstants)
        # λv, λe, λb = excitationfactor(ea, dFdθ, R, Rg, efconstants, rxcomponent)
        xmtrterm = λv*fzₜ*Cγ + λb*fyₜ*Sγ*Sϕ + λe*fxₜ*Sγ*Cϕ

        # Receiver term
        # TODO: Handle multiple components
        fzᵣ, fxᵣ, fyᵣ = heightgains(zᵣ, ea, frequency, efconstants)
        if rxcomponent == FC_Ez
            rcvrterm = fzᵣ
        elseif rxcomponent == FC_Ex
            rcvrterm = fxᵣ
        elseif rxcomponent == FC_Ey
            rcvrterm = fyᵣ
        end

        modesum += xmtrterm*rcvrterm*exp(-im*k*(S₀-1)*x)
    end

    Q = 682.2408*sqrt(frequency.f/1000*txpower/1000)  # in lw_sum_modes.for
    # Q = Z₀/(4π)*sqrt(2π*txpower/10k)*k/2  # Ferguson and Morfitt 1981 eq (21), V/m, NOT uV/m!
    # Q *= 100 # for V/m to uV/m

    # TODO: Radiation resistance correction if zₜ > 0

    E = Q/sqrt(abs(sin(x/Rₑ)))*modesum
    phase = angle(E)  # ranges between -π:π rad
    amp = 10log10(abs2(E))  # == 20log10(abs(E))

    return E, phase, amp
end

function Efield(
    modes,
    modeparams::ModeParameters,
    tx::Exciter,
    rx::GroundSampler
    )

    X = distance(rx)
    Xlength = length(X)

    Etype = eltype(modes[1])
    E = Array{Etype}(undef, Xlength)
    phase = Array{real(Etype)}(undef, Xlength)
    amp = Array{real(Etype)}(undef, Xlength)

    # TODO: This is really inefficient, because some of these things don't need to be computed
    # for each `x`, but for each `ea`, which is a much smaller size
    for i in eachindex(X)
        e, p, a = Efield(X[i], modes, modeparams, tx, rx)
        E[i] = e
        phase[i] = p
        amp[i] = a
    end

    unwrap!(phase)

    return E, phase, amp
end

function unwrap!(phasearray::AbstractArray)
    @inbounds for i in 2:length(phasearray)
        d = phasearray[i] - phasearray[i-1]
        if d >= π
            d -= 2π
        elseif d < -π
            d += 2π
        end
        phasearray[i] = phasearray[i-1] + d
    end
    return nothing
end

function referencetoground(ea::EigenAngle)
    return EigenAngle(asin(ea.sinθ/sqrt(1 - 2/Rₑ*H)))
end
