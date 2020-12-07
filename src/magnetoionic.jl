@doc raw"""
    susceptibility(altitude, frequency, bfield, species; params=LMPParams())

Compute the ionosphere susceptibility tensor `M` as a `SMatrix{3,3}` using
`species.numberdensity` and `species.collisionfrequency` at `altitude`.

If `params.earthcurvature == true`, `M` includes a first order correction for earth
curvature by means of a fictitious refractive index [^Pappert1967].

The susceptibility matrix is calculated from the constitutive relations presented in
[^Ratcliffe1959]. This includes the effect of earth's magnetic field vector and collisional
damping on electron motion.

The tensor is:
```math
M = -\frac{X}{U(U²-Y²)}
        \begin{pmatrix}
            U² - x²Y²    & -izUY - xyY² & iyUY - xzY² \\
            izUY - xyY²  & U² - y²Y²    & -ixUY - yzY² \\
            -iyUY - xzY² & ixUY - yzY²  & U² - z²Y²
        \end{pmatrix}
```
where ``X = ωₚ²/ω²``, ``Y = |ωₕ/ω|``, ``Z = ν/ω``, and ``U = 1 - iZ``. The earth curvature
correction subtracts ``2/Rₑ*(H - altitude)`` from the diagonal of ``M`` where ``H`` is
`params.curvatureheight`.

# References

[^Pappert1967]: R. A. Pappert, E. E. Gossard, and I. J. Rothmuller, “A numerical
    investigation of classical approximations used in VLF propagation,” Radio Science,
    vol. 2, no. 4, pp. 387–400, Apr. 1967, doi: 10.1002/rds196724387.

[^Ratcliffe1959]: J. A. Ratcliffe, "The magneto-ionic theory & its applications to the
    ionosphere," Cambridge University Press, 1959.
"""
function susceptibility(altitude, frequency, bfield, species; params=LMPParams())
    @unpack earthradius, earthcurvature, curvatureheight = params

    B, x, y, z = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    x², y², z² = x^2, y^2, z^2
    ω = frequency.ω

    e, m = species.charge, species.mass
    N, nu = species.numberdensity, species.collisionfrequency

    invω = inv(ω)
    invmω = invω/m  # == inv(m*ω)

    # Constitutive relations
    # (see Budden1955a, pg. 517, Budden1988 pg. 39, or [^Ratcliffe1959])
    X = N(altitude)*e^2*invω*invmω/E0
    Y = e*B*invmω  # [^Ratcliffe1959] pg. 182 specifies that sign of `e` is included here
    Z = nu(altitude)*invω
    U = 1 - 1im*Z

    # TODO: Support multiple species, e.g.
    # U²D = zero()
    # Y²D = zero()
    # UYD = zero()
    # for i = 1:length(species)
    #     U²D += U^2*D
    #     UYD += Y*U*D
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
    UYD = U*Y*D

    # Leverage partial symmetry of M to reduce computations
    izUYD = 1im*z*UYD
    xyY²D = x*y*Y²D
    iyUYD = 1im*y*UYD
    xzY²D = x*z*Y²D
    ixUYD = 1im*x*UYD
    yzY²D = y*z*Y²D

    # Elements of `M`
    M11 = U²D - x²*Y²D
    M21 = izUYD - xyY²D
    M31 = -iyUYD - xzY²D
    M12 = -izUYD - xyY²D
    M22 = U²D - y²*Y²D
    M32 = ixUYD - yzY²D
    M13 = iyUYD - xzY²D
    M23 = -ixUYD - yzY²D
    M33 = U²D - z²*Y²D

    if earthcurvature
        curvaturecorrection = 2/earthradius*(curvatureheight - altitude)

        M11 -= curvaturecorrection
        M22 -= curvaturecorrection
        M33 -= curvaturecorrection
    end

    M = SMatrix{3,3}(M11, M21, M31, M12, M22, M32, M13, M23, M33)

    return M
end
"""
    susceptibility(altitude, me::ModeEquation; params=LMPParams())
"""
susceptibility(altitude, me::ModeEquation; params=LMPParams()) =
    susceptibility(altitude, me.frequency, me.waveguide, params=params)
"""
    susceptibility(altitude, frequency, w::HomogeneousWaveguide; params=LMPParams())
"""
susceptibility(altitude, frequency, w::HomogeneousWaveguide; params=LMPParams()) =
    susceptibility(altitude, frequency, w.bfield, w.species, params=params)


"""
    tmatrix(ea::EigenAngle, M)

Return the matrix components of `T` whose elements are used in the calculation
of matrix `S` for the differential equations for the ionospheric reflection coefficient `R`.

Following Budden's formalism [^Budden1955a] the differential equations for calculating the
reflection coefficient for a radio wave obliquely incident on a horizontally stratified
ionosphere are derived from Maxwell's equations in conjunction with the constitutive relations
for the ionosphere, represented by the matrix `M`. After eliminating the vertically directed
components ``Ez`` and ``Hz``, the equations can be written ``∂e/∂s = -i T e`` where ``s`` is
a proxy for height `z`, and ``e = (Ex -Ey Hx Hy)``, as detailed by [^Clemmow1954]. This
function calculates components of the matrix `T`.

See also: [`susceptibility`](@ref)

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing
    reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227,
    no. 1171, pp. 516–537, Feb. 1955.

[^Clemmow1954]: P. C. Clemmow and J. Heading, “Coupled forms of the differential equations
    governing radio propagation in the ionosphere,” Mathematical Proceedings of the
    Cambridge Philosophical Society, vol. 50, no. 2, pp. 319–333, Apr. 1954.
"""
function tmatrix(ea::EigenAngle, M)
    S, C² = ea.sinθ, ea.cos²θ

    # Denominator of most of the entries of `T`
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

    # `TMatrix` is a special 4×4 matrix
    return TMatrix(T11, T31, T41,
                   T12, T32, T42,
                   T14, T34, T44)
end

"""
    tmatrix(ea::EigenAngle, M, ::Dθ)

Return a dense `SMatrix` with the derivative of `T` with respect to `θ` at eigenangle `ea`.
"""
function tmatrix(ea::EigenAngle, M, ::Dθ)
    S, C = ea.sinθ, ea.cosθ
    dC² = -2*S*C  # d/dθ (C²)

    den = inv(1 + M[3,3])

    dT11 = -C*M[3,1]*den
    dT12 = C*M[3,2]*den
    dT13 = 0
    dT14 = dC²*den
    dT21 = 0
    dT22 = 0
    dT23 = 0
    dT24 = 0
    dT31 = 0
    dT32 = dC²
    dT33 = 0
    dT34 = C*M[2,3]*den
    dT41 = 0
    dT42 = 0
    dT43 = 0
    dT44 = -C*M[1,3]*den

    return SMatrix{4,4}(dT11, dT21, dT31, dT41,
                        dT12, dT22, dT32, dT42,
                        dT13, dT23, dT33, dT43,
                        dT14, dT24, dT34, dT44)
end

"""
    bookerquartic(ea, M)

Compute roots `q` and the coefficients `B` of the Booker quartic for `EigenAngle` `ea` and
susceptibility tensor `M`.

This function uses `PolynomialRoots.jl` to find the roots.

See also: [`dbookerquartic`](@ref)

# References

[^Sheddy1968a]: C. H. Sheddy, “A General Analytic Solution for Reflection From a Sharply
    Bounded Anisotropic Ionosphere,” Radio Science, vol. 3, no. 8, pp. 792–795, Aug. 1968.
"""
function bookerquartic(ea::EigenAngle, M)
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
    B2 = -C²pM33*M11p1 + M13M31 - M33p1*C²pM22 + M23M32
    B1 = S*(M12M23 + M21M32 - C²pM22*M13pM31)
    B0 = M11p1*C²pM22*C²pM33 +
         M12M23*M[3,1] + M[1,3]*M21M32 -
         M13M31*C²pM22 - M11p1*M23M32 -
         M[1,2]*M[2,1]*C²pM33

    booker_quartic_coeffs = MVector{5}(B0, B1, B2, B3, B4)
    booker_quartic_roots = MVector{4,eltype(booker_quartic_coeffs)}(0, 0, 0, 0)

    # `roots!` dominates this functions runtime
    roots!(booker_quartic_roots, booker_quartic_coeffs, NaN, 4, false)

    return booker_quartic_roots, SVector(booker_quartic_coeffs)
end

"""
    bookerquartic(T::TMatrix)

Solve the Booker quartic in depressed form in terms of `T`.
"""
function bookerquartic(T::TMatrix)
    # This technique is used in LWPC's wavefield subroutines, e.g. "wf_init.for".

    # Precompute
    T34T42 = T[3,4]*T[4,2]

    # This is the depressed form of the quartic
    b3 = -T[1,1] - T[4,4]
    b2 = T[1,1]*T[4,4] - T[1,4]*T[4,1] - T[3,2]
    b1 = T[3,2]*(T[1,1] + T[4,4]) - T[1,2]*T[3,1] - T34T42
    b0 = -T[1,1]*(T[3,2]*T[4,4] - T34T42) +
        T[1,2]*(T[3,1]*T[4,4] - T[3,4]*T[4,1]) -
        T[1,4]*(T[3,1]*T[4,2] - T[3,2]*T[4,1])

    booker_quartic_coeffs = MVector{5}(b0, b1, b2, b3, 1)
    booker_quartic_roots = MVector{4,eltype(booker_quartic_coeffs)}(0, 0, 0, 0)

    roots!(booker_quartic_roots, booker_quartic_coeffs, NaN, 4, false)

    return booker_quartic_roots, SVector(booker_quartic_coeffs)
end

"""
    dbookerquartic(ea::EigenAngle, M, q, B)

Compute derivative `dq` of the Booker quartic roots with respect to ``θ``.

See also: [`bookerquartic`](@ref)
"""
function dbookerquartic(ea::EigenAngle, M, q, B)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    dS = C
    dC = -S
    dC² = -2*S*C

    dB3 = dS*(M[1,3] + M[3,1])
    dB2 = -dC²*(2 + M[1,1] + M[3,3])
    dB1 = dS/S*B[2] - S*dC²*(M[1,3] + M[3,1])
    dB0 = dC²*(2*C²*(1 + M[1,1]) + M[3,3] + M[2,2] + M[1,1]*(M[3,3] + M[2,2]) -
            M[1,3]*M[3,1] - M[1,2]*M[2,1])

    dq = similar(q)
    @inbounds for i in eachindex(dq)
        dq[i] = -(((dB3*q[i] + dB2)*q[i] + dB1)*q[i] + dB0) /
                (((4*B[5]*q[i] + 3*B[4])*q[i] + 2*B[3])*q[i] + B[2])
    end

    return SVector(dq)
end

"""
    dbookerquartic(T::TMatrix, dT, q, B)
"""
function dbookerquartic(T::TMatrix, dT, q, B)
    dB3 = -dT[1,1] - dT[4,4]
    dB2 = T[1,1]*dT[4,4] + dT[1,1]*T[4,4] - dT[1,4]*T[4,1] - dT[3,2]
    dB1 = -(-T[3,2]*(dT[1,1] + dT[4,4]) - dT[3,2]*(T[1,1] + T[4,4]) + dT[1,2]*T[3,1] +
        dT[3,4]*T[4,2])
    dB0 = -T[1,1]*(T[3,2]*dT[4,4] + dT[3,2]*T[4,4] - dT[3,4]*T[4,2]) -
        dT[1,1]*(T[3,2]*T[4,4] - T[3,4]*T[4,2]) +
        T[1,2]*(T[3,1]*dT[4,4] - dT[3,4]*T[4,1]) + dT[1,2]*(T[3,1]*T[4,4] - T[3,4]*T[4,1]) -
        T[1,4]*(-dT[3,2]*T[4,1]) - dT[1,4]*(T[3,1]*T[4,2] - T[3,2]*T[4,1])

    dq = similar(q)
    @inbounds for i in eachindex(dq)
        dq[i] = -(((dB3*q[i] + dB2)*q[i] + dB1)*q[i] + dB0) /
                (((4*B[5]*q[i] + 3*B[4])*q[i] + 2*B[3])*q[i] + B[2])
    end

    return SVector(dq)
end

"""
    upgoing(v)

Calculate the absolute angle of `v` in radians from 315°×π/180 on the complex plane. Smaller
values indicate upgoing waves.
"""
function upgoing(v::Complex)
    # -π/4 is 315°
    a = -π/4 - angle(v)

    # angle() is +/-π
    # `a` can never be > 3π/4, but there is a region on the plane where `a` can be < -π
    a < -π && (a += 2π)
    return abs(a)
end

upgoing(v::Real) = oftype(v, π/4)

"""
    sortquarticroots!(q)

Sort array of quartic roots `q` in place such that the first two correspond to upgoing waves and the
latter two correspond to downgoing waves.

General locations of the four roots of the Booker quartic on the complex plane
corresponding to the:

    1) upgoing evanescent wave
    2) upgoing travelling wave
    3) downgoing evanescent wave
    4) downgoing travelling wave

```
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
```

Based on [^Pitteway1965] fig. 5.

[^Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields, reflexion
    coefficients and polarizations for long radio waves in the lower ionosphere. I.,” Phil.
    Trans. R. Soc. Lond. A, vol. 257, no. 1079, pp. 219–241, Mar. 1965,
    doi: 10.1098/rsta.1965.0004.
"""
function sortquarticroots!(q)
    # It looks like we could just `sort!(q, by=real)`, but the quadrants of each root are
    # not fixed and the positions in the Argand diagram are just approximate.

    length(q) != 4 && @warn "length of `q` is not 4"

    # Calculate and sort by distance from 315°
    # The two closest are upgoing and the two others are downgoing
    # This is faster than `sort!(q, by=upgoing)`
    dist = MVector{4,real(eltype(q))}(undef)
    @inbounds for i in eachindex(q)
        dist[i] = upgoing(q[i])
    end

    i = 4
    @inbounds while i > 1
        if dist[i] < dist[i-1]
            dist[i], dist[i-1] = dist[i-1], dist[i]
            q[i], q[i-1] = q[i-1], q[i]
            i = 4  # Need to restart at the end
        else
            i -= 1
        end
    end

    # Now arrange the roots so that they correspond to:
    # 1) upgoing evanescent, 2) upgoing travelling
    # 3) downgoing evanescent, and 4) downgoing travelling
    if imag(q[1]) > imag(q[2])
        q[2], q[1] = q[1], q[2]
    end
    if imag(q[3]) < imag(q[4])
        q[4], q[3] = q[3], q[4]
    end

    # NOTE: Nagano et al 1975 says that of q1 and q2, the one with the largest
    # absolute value corresponds to e1 and the other to e2, although this method
    # does not guarantee that. In fact, I find that this criteria is often not met.

    return q
end
