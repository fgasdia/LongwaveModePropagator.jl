
"""
    bookerquartic(ea, M)

Return roots `q` and the coefficients `B` of the Booker quartic.

The Booker quartic is used in the solution of `R` for a sharply bounded ionosphere. This
function uses the `PolynomialRoots` package to find the roots.

See also: [`sharplybounded_R`](@ref)

# References

[^Sheddy1968a]: C. H. Sheddy, “A General Analytic Solution for Reflection From a Sharply Bounded Anisotropic Ionosphere,” Radio Science, vol. 3, no. 8, pp. 792–795, Aug. 1968.
"""
function bookerquartic!(ea::EigenAngle, M)
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
function bookerquartic!(T::TMatrix)
    # This is the depressed form of the quartic
    b3 = -(T[1,1] + T[4,4])
    b2 = T[1,1]*T[4,4] - T[1,4]*T[4,1] - T[3,2]
    b1 = -(-T[3,2]*(T[1,1] + T[4,4]) + T[1,2]*T[3,1] + T[3,4]*T[4,2])
    b0 = -T[1,1]*(T[3,2]*T[4,4] - T[3,4]*T[4,2]) +
        T[1,2]*(T[3,1]*T[4,4] - T[3,4]*T[4,1]) -
        T[1,4]*(T[3,1]*T[4,2] - T[3,2]*T[4,1])

    BOOKER_QUARTIC_COEFFS.data = (b0, b1, b2, b3, complex(1.0))

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

    # NOTE: Nagano et al 1975 says that of v1 and v2, the one with the largest
    # absolute value corresponds to e1 and the other to e2, although this method
    # does not guarantee that. In fact, I find that this criteria is often not met.

    return v
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
function susceptibility(z, frequency::Frequency, bfield::BField, species::Species)
    B, lx, ly, lz = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    lx², ly², lz² = lx^2, ly^2, lz^2
    ω = frequency.ω

    # Constitutive relations (see Budden1955a, pg. 517 or Budden1988 pg. 39)
    e, m, N, ν = species.charge, species.mass, species.numberdensity, species.collisionfrequency
    invω = inv(ω)
    invmω = invω/m  # == inv(m*ω)

    X = N(z)*e^2*invω*invmω/E0  # ωₚ²/ω² plasma frequency / ω
    Y = abs(e*B*invmω)  # |ωₕ/ω|  gyrofrequency / ω  # Nagano et al 1975 specifies |ωₕ/ω|
    Z = ν(z)*invω  # collision frequency / ω
    U = complex(1, -Z)  # == 1 - im*Z

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

    earthcurvature = 2/EARTHRADIUS*(CURVATURE_HEIGHT - z)

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
    M = SMatrix{3,3,ComplexF64,9}(M11, M21, M31,
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
