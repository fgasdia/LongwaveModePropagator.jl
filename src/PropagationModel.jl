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

# Not great, but can be changed as `EARTHCURVATURE[]=false`
# TODO: where does this need to be considered?
const EARTHCURVATURE = Ref(true)

# Passing MArrays between functions causes allocations. They are avoided by
# mutating this const in place. `roots!` requires Complex values.
const BOOKER_QUARTIC_ROOTS = MVector{4}(zeros(ComplexF64, 4))
const BOOKER_QUARTIC_COEFFS = MVector{5,ComplexF64}(undef)

# const BOOKER_QUARTIC_ROOTS = MVector{4}(zeros(Complex{BigFloat}, 4))
# const BOOKER_QUARTIC_COEFFS = MVector{5,Complex{BigFloat}}(undef)

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
    Rg::SDiagonal{2,T}
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

"""
    WavefieldIntegrationParams{T1,T2,T3,F,G}(z, ortho_scalar, e1_scalar, e2_scalar, ea, frequency, bfield, species)

Parameters passed to Pitteway integration of wavefields.

    `z::T1`
    `ortho_scalar::Complex{T2}`
    `e1_scalar::T2`
    `e2_scalar::T2`
    `ea::EigenAngle{T3}`
    `frequency::Frequency`
    `bfield::BField`
    `species::Constituent{F,G}`
"""
@with_kw struct WavefieldIntegrationParams{T1,T2,T3,F,G}
    z::T1
    ortho_scalar::Complex{T2}
    e1_scalar::T2
    e2_scalar::T2
    ea::EigenAngle{T3}
    frequency::Frequency
    bfield::BField
    species::Constituent{F,G}
end

"""
    WavefieldIntegrationParams{T}(ea, frequency, bfield, species)

Initialize a `WavefieldIntegrationParams` for downward Pitteway scaled
integration. Requires the parameter `T`, which should be the wavefield type
(usually `Complex{Float64}`).

Automatically set values are:

    `z = TOPHEIGHT`
    `ortho_scalar = zero(complex(T2))`
    `e1_scalar = one(real(T2))`
    `e2_scalar = one(real(T2))`

`ortho_scalar` will always be complex and `e1_scalar` and `e2_scalar` will
always be real, so it is sufficient for `Float64` to be provided as `T` even
for complex wavefields.
"""
function WavefieldIntegrationParams{T}(ea::EigenAngle{T3}, frequency::Frequency, bfield::BField, species::Constituent{F,G}) where {T,T3,F,G}
    return WavefieldIntegrationParams{typeof(TOPHEIGHT),real(T),T3,F,G}(TOPHEIGHT, zero(complex(T)), one(real(T)), one(real(T)), ea, frequency, bfield, species)
end

"""
    ScaleRecord{T1,T2}(z, e, ortho_scalar, e1_scalar, e2_scalar)

Struct used for saving wavefield scaling information during Pitteway integration
of wavefields.

!!! note

    `T2` must be real, so `ScaleRecord` would have type, e.g.
    `ScaleRecord{eltype(zs), real(eltype(e0))}`
"""
struct ScaleRecord{T1,T2<:Real}
    z::T1
    e::SMatrix{4,2,Complex{T2},8}
    ortho_scalar::Complex{T2}
    e1_scalar::T2
    e2_scalar::T2
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
function bookerquartic!(ea::EigenAngle, M::AbstractArray{eltype(BOOKER_QUARTIC_ROOTS)})
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
function bookerquartic!(T::TMatrix{eltype(BOOKER_QUARTIC_COEFFS)})
    # This is the depressed form of the quartic
    b3 = -(T[1,1] + T[4,4])
    b2 = T[1,1]*T[4,4] - T[1,4]*T[4,1] - T[3,2]
    b1 = -(-T[3,2]*(T[1,1] + T[4,4]) + T[1,2]*T[3,1] + T[3,4]*T[4,2])
    b0 = -T[1,1]*(T[3,2]*T[4,4] - T[3,4]*T[4,2]) +
        T[1,2]*(T[3,1]*T[4,4] - T[3,4]*T[4,1]) -
        T[1,4]*(T[3,1]*T[4,2] - T[3,2]*T[4,1])

    BOOKER_QUARTIC_COEFFS.data = (b0, b1, b2, b3, one(eltype(BOOKER_QUARTIC_COEFFS)))

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

General locations of the four roots of the Booker quartic on the complex plane
corresponding to the:

    1) upgoing evanescent wave
    2) upgoing travelling wave
    3) downgoing evanescent wave
    4) downgoing travelling wave

From Pitteway 1967 fig 5. Sheddy 1968 also mentions that the two upgoing waves
correspond to the two solutions that lie close to the positive real axis and the
negative imaginary axis, although in [`sharpboundaryreflection`](@ref) the order
of the two upgoing waves doesn't matter---swapping q₁ and q₂ will result in the
same reflection coefficient matrix.

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

    # Uncomment to test if order of q₁, q₂ matters for `sharpboundaryreflection()`
    # v[2], v[1] = v[1], v[2]

    # BUG: (?) Nagano et al 1975 says that of v1 and v2, the one with the largest
    # absolute value corresponds to e1 and the other to e2, although this method
    # does not guarantee that. In fact, I find that this criteria is often not met.

    return v
end

function _sharpboundaryreflection(ea::EigenAngle, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # XXX: `bookerquartic` (really `roots!`) dominates this functions runtime
    bookerquartic!(ea, M)

    # We choose the 2 roots corresponding to upward travelling waves as being
    # those that lie close to the positive real axis and negative imaginary axis
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

    # NOTE: `BOOKER_QUARTIC_ROOTS` is sorted in _sharpboundaryreflection
    # Could write an `issorted` function, but sorting time is dominated by
    # `upgoing`, which would presumably be required by an `issorted`
    q = BOOKER_QUARTIC_ROOTS

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
    Y = abs(e*B/mω)  # |ωₕ/ω|  gyrofrequency / ω  # Nagano et al 1975 specifies |ωₕ/ω|
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
    M11 = U²D - lx²*Y²D
    M21 = izYUD - xyY²D
    M31 = -iyYUD - xzY²D
    M12 = -izYUD - xyY²D
    M22 = U²D - ly²*Y²D
    M32 = ixYUD - yzY²D
    M13 = iyYUD - xzY²D
    M23 = -ixYUD - yzY²D
    M33 = U²D - lz²*Y²D

    if EARTHCURVATURE[]
        M11 -= earthcurvature
        M22 -= earthcurvature
        M33 -= earthcurvature
    end

    # Remember, column major
    M = SMatrix{3,3}(M11, M21, M31,
                     M12, M22, M32,
                     M13, M23, M33)

    return M
end

"""
    tmatrix(ea, M)

Return the matrix components of `T` whose elements are used in the calculation
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

    # Remember, column major. And it's actually a special 4×4.
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

    # TODO: Integration method
    # TODO: Does saving actually alter performance?
    prob = ODEProblem{false}(dRdz, Rtop, (TOPHEIGHT, BOTTOMHEIGHT), (ea, modeparams))

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8, save_on=false, save_end=true)

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

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8, save_on=false, save_end=true)

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

    Rg = SDiagonal(Rg11, Rg22)

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
function findmodes(origcoords::AbstractArray{T}, modeparams::ModeParameters, tolerance=1e-6) where {T}
    zroots, zpoles = grpf(θ->solvemodalequation(EigenAngle{T}(θ), modeparams),
                          origcoords, tolerance, 30000)

    return zroots
    # return @SVector [EigenAngle{T}(r) for r in zroots]
end

################
# Wavefields
################

"""
    initialwavefields(T::TMatrix)

Calculate the initial wavefields vector ``[Ex₁ Ex₂
                                           Ey₁ Ey₂
                                           ℋx₁ ℋx₂
                                           ℋy₁ ℋy₂]``
for the two upgoing wavefields where subscript `1` is the evanescent wave and
`2` is the travelling wave.

This function solves the equation ``Te = qe``, equivalently the eigenvector
problem ``(T- qI)e = 0``. First, the Booker Quartic is solved for the roots `q`,
and they are sorted so that the roots associated with the two upgoing waves are
selected, where eigenvalue ``q₁`` corresponds to the evanescent wave and ``q₂``
the travelling wave. Then `e` is solved as the eigenvectors for the two `q`s. An
analytical solution is used where `e[2,:] = 1`.
"""
function initialwavefields(T::TMatrix)
    # TODO: rename to bookerwavefields?

    bookerquartic!(T)
    sortquarticroots!(BOOKER_QUARTIC_ROOTS)

    q = BOOKER_QUARTIC_ROOTS

    # Precompute
    T14T41 = T[1,4]*T[4,1]
    T14T42 = T[1,4]*T[4,2]
    T12T41 = T[1,2]*T[4,1]

    # Temporary MArray for filling in wavefields
    # (04/2020) Somehow, this is slightly faster than hard coding the 1s and 2s
    e = MArray{Tuple{4,2},eltype(BOOKER_QUARTIC_ROOTS)}(undef)

    for i = 1:2
        d = T14T41 - (T[1,1] - q[i])*(T[4,4] - q[i])
        dinv = 1/d

        e[1,i] = (T[1,2]*(T[4,4] - q[i]) - T14T42)*dinv
        e[2,i] = 1
        e[3,i] = q[i]
        e[4,i] = (-T12T41 + T[4,2]*(T[1,1] - q[i]))*dinv
    end

    # By returning as SArray instead of MArray, the MArray doesn't get hit by GC
    return SArray(e)
end

"""
    dedz(e, k, T::Tmatrix)

Calculates derivative of field components vector ``de/dz = -i k T e``.
"""
dedz(e, k, T::TMatrix) = -im*k*(T*e)  # `(T*e)` uses specialized TMatrix math

"""
    dedz(e, p, z)

Calculates derivative of field components vector `e` at height `z`.

The parameters tuple `p` should contain (`Frequency`, `BField`, `Constituent`)
or be a `WavefieldIntegrationParams` struct. This function internally calls
[`susceptibility`](@ref) and [`tmatrix`](@ref) and is typically used by
[`integratewavefields`](@ref).
"""
function dedz(e, p, z)
    @unpack ea, frequency, bfield, species = p

    M = susceptibility(z, frequency, bfield, species)
    T = tmatrix(ea, M)

    return dedz(e, frequency.k, T)
end

"""
    scalingcondition(e, z, integrator)

Return `true` if wavefields should be scaled, otherwise `false`.

Specifically, if any component of `real(e)` or `imag(e)` are `>= 1`, return
`true`.
"""
scalingcondition(e, z, integrator) = any(x -> (real(x) >= 1 || imag(x) >= 1), e)

"""
    scale!(integrator)

Apply wavefield scaling with [`scalewavefields`](@ref) to the integrator.
"""
function scale!(integrator)
    new_e, new_orthos, new_e1s, new_e2s = scalewavefields(integrator.u)

    # Last set of scaling values
    @unpack ea, frequency, bfield, species = integrator.p

    #==
    NOTE: `integrator.t` is the "time" of the _proposed_ step. Therefore,
    integrator.t` might equal `0.0`, for example, before it's actually gotten
    to the bottom. `integrator.prevt` is the last `t` on the "left"
    side of the `integrator`, which covers the local interval [`tprev`, `t`].
    The "condition" is met at `integrator.t` and `integrator.t` is the time at
    which the affect occurs.
    However, it is not guaranteed that each (`tprev`, `t`) directly abuts the
    next `tprev`, `t`).
    ==#

    # NOTE: we must entirely reconstruct the entire NamedTuple from scratch
    integrator.p = WavefieldIntegrationParams(integrator.t,
                                              new_orthos,
                                              new_e1s, new_e2s,
                                              ea, frequency, bfield, species)

    integrator.u = new_e

    return nothing
end

"""
    save_values(u, t, integrator)

Return a `ScaleRecord` from `u`, `t`, and `integrator`.

Used by SavingCallback in [`integratewavefields`](@ref).
"""
save_values(u, t, integrator) = ScaleRecord(integrator.p.z,
                                            u,
                                            integrator.p.ortho_scalar,
                                            integrator.p.e1_scalar,
                                            integrator.p.e2_scalar)

"""
    scalewavefields(e1, e2)

Returns orthonormalized vectors `e1` and `e2`, as well as the scaling terms `a`,
`e1_scale_val`, and `e2_scale_val` applied to the original vectors.

First applies Gram-Schmidt orthogonalization and then scales the vectors so they
each have length 1, i.e. `norm(e1) == norm(e2) == 1`. This is the technique
suggested by [^Pitteway1965] to counter numerical swamping during integration of
wavefields.

# References

[^Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields,
reflexion coefficients and polarizations for long radio waves in the lower
ionosphere. I.,” Phil. Trans. R. Soc. Lond. A, vol. 257, no. 1079,
pp. 219–241, Mar. 1965, doi: 10.1098/rsta.1965.0004.
"""
function scalewavefields(e1::AbstractVector, e2::AbstractVector)
    # Orthogonalize vectors `e1` and `e2` (Gram-Schmidt process)
    # `dot` for complex vectors automatically conjugates first vector
    e1_dot_e1 = real(dot(e1, e1))  # == sum(abs2.(e1)), `imag(dot(e1,e1)) == 0`
    a = dot(e1, e2)/e1_dot_e1  # purposefully unsigned XXX: necessary?
    e2 -= a*e1

    # Normalize `e1` and `e2`
    e1_scale_val = 1/sqrt(e1_dot_e1)
    e2_scale_val = 1/norm(e2)  # == 1/sqrt(dot(e2,e2))
    e1 *= e1_scale_val
    e2 *= e2_scale_val  # == normalize(e2)

    return e1, e2, a, e1_scale_val, e2_scale_val
end

"""
    scalewavefields(e)

!!! note

    This function only applies scaling to the first 2 columns of `e`.
"""
function scalewavefields(e::AbstractArray)
    e1, e2, a, e1_scale_val, e2_scale_val = scalewavefields(e[:,1], e[:,2])

    if size(e, 2) == 2
        e = hcat(e1, e2)
    else
        e = hcat(e1, e2, e[:,3:end])
    end

    return e, a, e1_scale_val, e2_scale_val
end

"""
    unscalewavefields!(e, saved_values::SavedValues)

Unscale the integrated wavefields `e` in place.

Assumes fields have been scaled by [`scalewavefields`](@ref) during integration.

See also [`unscalewavefields`](@ref) for additional details.
"""
function unscalewavefields!(e::AbstractVector, saved_values::SavedValues)

    zs = saved_values.t
    records = saved_values.saveval

    # Initialize the "reference" scaling altitude
    ref_z = last(records).z

    # Usually `osum = 0`, `prod_e1 = 1`, and `prod_e2 = 1` at initialization,
    # but we set up the fields at the ground (`last(records)`) outside the loop
    osum = last(records).ortho_scalar
    prod_e1 = last(records).e1_scalar
    prod_e2 = last(records).e2_scalar

    # Unscaling we go from the bottom up
    @inbounds for i in reverse(eachindex(e))

        # Unpack variables
        record_z = records[i].z
        scaled_e = records[i].e
        ortho_scalar = records[i].ortho_scalar
        e1_scalar = records[i].e1_scalar
        e2_scalar = records[i].e2_scalar

        z = zs[i]

        # TODO: `zs[i]` and `records[i].e` dominate the entire funtcion runtime
        # Due to iterating in reverse? Need to investigate. Possibly reverse the
        # arrays first but would then need to reverse output

        # Only update ref_z when there is both:
        # 1) we have reached the height where a new ref_z should go into effect
        # 2) there is a new ref_z
        if (z > record_z) & (record_z > ref_z)
            ref_z = record_z

            osum *= e1_scalar/e2_scalar
            osum += ortho_scalar
            prod_e1 *= e1_scalar
            prod_e2 *= e2_scalar
        end

        if i == lastindex(e)  # == [end]
            # Bottom doesn't require correction
            e[i] = scaled_e
        elseif z > ref_z
            # From the bottom, the first correction may not need to be applied
            # until some higher altitude
            e2 = (scaled_e[:,2] - osum*scaled_e[:,1])*prod_e2
            e1 = scaled_e[:,1]*prod_e1
            e[i] = hcat(e1,e2)
        else
            e[i] = scaled_e
        end
    end

    return nothing
end

"""
    unscalewavefields(saved_values::SavedValues)

Return the unscaled integrated wavefields originally scaled by
[`scalewavefields`](@ref).

The fields are not "unscaled" so much as further scaled so that at all heights
the fields are scaled similarly to the cumulative scaling that applied to the
fields at the bottom.

The bottom level does not get unscaled. We reference the higher levels to the
bottom. The level above the bottom level needs to be additionally scaled by the
amount that was applied to originally get from this level down to the bottom
level. The next level up (2 above the bottom level) needs to be scaled by the
amount applied to the next level and then the bottom level, i.e. we keep track
of a cumulative correction on the way back up.
"""
function unscalewavefields(saved_values::SavedValues)
    # Array of SArray{Tuple{4,2}, Complex{Float64}}
    e = Vector{typeof(saved_values.saveval[1].e)}(undef, length(saved_values.saveval))

    unscalewavefields!(e, saved_values)

    return e
end

"""

Always integrates to ground (0).

zs should be going down (although this isn't strictly enforced by this function)
"""
function integratewavefields(zs, ea, frequency, bfield, species)
    # TODO: version that updates output `e` in place

    # Initial conditions
    Mtop = susceptibility(first(zs), frequency, bfield, species)
    Ttop = tmatrix(ea, Mtop)
    e0 = initialwavefields(Ttop)

    # Normalize e0 (otherwise top fields are out of whack with unscaled fields
    # because scaling at top is unity).
    # Don't want to orthogonalize `e2` here because it won't be undone!
    e0 = hcat(normalize(e0[:,1]), normalize(e0[:,2]))

    # saved_positions=(true, true) because we discontinuously modify `u`. This is
    # independent of saveat and save_everystep
    cb = DiscreteCallback(scalingcondition, scale!, save_positions=(true, true))

    saved_values = SavedValues(eltype(zs), ScaleRecord{eltype(zs), real(eltype(e0))})

    # `save_everystep` because we need to make sure we save when affect! occurs
    # `saveat=zs[2:end-1]` because otherwise we double save end points
    scb = SavingCallback(save_values, saved_values,
                         save_everystep=true, saveat=zs[2:end-1],
                         tdir=sign(last(zs)-first(zs)))


    p = WavefieldIntegrationParams{eltype(e0)}(ea, frequency, bfield, species)

    # (May 5, 2020) DifferentialEquations chooses Vern9(false) on daytime Wait
    # profile. But Vern6, 7, or 8 are likely more efficient since the stiff
    # Rodas5 algorithm was never called.

    # WARNING: Without `lazy=false` (since we're using DiscreteCallback) don't
    # use continuous solution output!
    prob = ODEProblem{false}(dedz, e0, (first(zs), last(zs)), p)
    sol = solve(prob, callback=CallbackSet(cb, scb),
                save_everystep=false, save_start=false, save_end=false)

    e = unscalewavefields(saved_values)

    # TODO: Only return at zs (not the extra points at integration steps)
    # then don't need to return saved_values.t
    return (e, saved_values.t)
end

integratewavefields(ea, frequency, bfield, species) = integratewavefields(TOPHEIGHT:-100:BOTTOMHEIGHT, ea, frequency, bfield, species)

"""
    vacuumreflectioncoeffs(ea, e)

Return ionosphere reflection coefficient matrix from upgoing wave fields `e`.

Integrating for one set of horizontal field components ``e = (Ex, -Ey, Z₀Hx, Z₀Hy)ᵀ``
can be separated into an upgoing and downgoing wave, each of which is generally
elliptically polarized. One might assume that the ratio of the amplitudes of
these two waves would give a reflection coefficient, and it does, except the
coefficient would only apply for an incident wave of that particular elliptical
polarization. However, the first set of fields can be linearly combined with
a second independent solution for the fields, which will generally have a
different elliptical polarization than the first. Two linear combinations of the
two sets of fields are formed with unit amplitude, linearly polarized
incident waves. The reflected waves then give the components ``R₁₁``, ``R₂₁`` or
``R₁₂``, ``R₂₂`` for the incident wave in the plane of incidence and
perpendicular to it, respectively [^Budden1988] (pg 552).

The process for determining the reflection coefficient requires resolving the
two sets of fields `e1` and `e2` into the four linearly polarized vacuum
modes. The layer of vacuum can be assumed to be so thin that it does not affect
the fields. There will be two upgoing waves, one of which has ``E``, and the
other ``H`` in the plane of incidence, and two downgoing waves, with ``E`` and
``H`` in the plane of incidence. If ``f₁, f₂, f₃, f₄`` are the complex
amplitudes of the four component waves, then in matrix notation ``e = Sᵥ f``.

For `e1` and `e2`, we can find the corresponding vectors `f1` and `f2` by
``f1 = Sᵥ⁻¹ e1``, ``f2 = Sᵥ⁻¹ e2`` where the two column vectors are partitioned
such that ``f1 = (u1, d1)ᵀ`` and ``f2 = (u2, d2)ᵀ`` for upgoing and downgoing
2-element vectors `u` and `d`. From the definition of the reflection coefficient
`R`, ``d = Ru``. Letting ``U = (u1, u2)``, ``D = (d1, d2)``, then ``D = RU`` and
the reflection coefficient is ``R = DU¹``. Because the reflection coefficient
matrix is a ratio of fields, either `e1` and/or `e2` can be independently
multiplied by an arbitrary constant and the value of `R` is unaffected.

For additional details, see [^Budden1988], chapter 18, section 7.

# References

[^Budden1988] K. G. Budden, The propagation of radio waves: the theory of radio
waves of low power in the ionosphere and magnetosphere, First paperback edition.
New York: Cambridge University Press, 1988.
"""
function vacuumreflectioncoeffs(ea::EigenAngle{T}, e1::AbstractArray{T2}, e2::AbstractArray{T2}) where {T,T2}
    C = ea.cosθ
    Cinv = 1/C

    # TODO: Special Sv matrix (also useful elsewhere?)
    Sv_inv = SMatrix{4,4,T,16}(Cinv, 0, -Cinv, 0,
                               0, -1, 0, -1,
                               0, -Cinv, 0, Cinv,
                               1, 0, 1, 0)

    f1 = Sv_inv*e1
    f2 = Sv_inv*e2

    out_type = promote_type(T, T2)
    U = SMatrix{2,2,out_type,4}(f1[1], f1[2], f2[1], f2[2])
    D = SMatrix{2,2,out_type,4}(f1[3], f1[4], f2[3], f2[4])

    return D/U
end

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
