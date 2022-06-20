#==
Functions related to the Booker quartic equation and the four waves associated with the
quartic roots.

Besides calculating the quartic roots themselves, this file also includes functions
related to calculating wavefields and vacuum reflection coefficients from the quartic roots.
==#

"""
    bookerquartic(ea, M)
    bookerquartic(T::TMatrix)

Compute roots `q` and the coefficients `B` of the Booker quartic described by the
susceptibility tensor `M` or `T` matrix.

# References

[Budden1988]: K. G. Budden, “The propagation of radio waves: the theory of radio
    waves of low power in the ionosphere and magnetosphere,” First paperback edition.
    New York: Cambridge University Press, 1988.
"""
bookerquartic

function bookerquartic(ea, M)
    S, C = sincos(ea)
    C² = C^2

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

function bookerquartic(T::TMatrix)
    # Precompute
    T34T42 = T[3,4]*T[4,2]

    # This is the depressed form of the quartic
    b3 = -T[1,1] - T[4,4]
    b2 = T[1,1]*T[4,4] - T[1,4]*T[4,1] - T[3,2]
    b1 = T[3,2]*(T[1,1] + T[4,4]) - T[1,2]*T[3,1] - T34T42
    b0 = T[1,1]*(T34T42 - T[3,2]*T[4,4]) +
        T[1,2]*(T[3,1]*T[4,4] - T[3,4]*T[4,1]) +
        T[1,4]*(T[3,2]*T[4,1] - T[3,1]*T[4,2])

    booker_quartic_coeffs = MVector{5}(b0, b1, b2, b3, 1)
    booker_quartic_roots = MVector{4,eltype(booker_quartic_coeffs)}(0, 0, 0, 0)

    roots!(booker_quartic_roots, booker_quartic_coeffs, NaN, 4, false)

    return booker_quartic_roots, SVector(booker_quartic_coeffs)
end

"""
    dbookerquartic(ea, M, q, B)
    dbookerquartic(T::TMatrix, dT, q, B)

Compute derivative `dq` of the Booker quartic roots `q` with respect to ``θ`` for the
ionosphere described by susceptibility tensor `M` or `T` matrix.
"""
dbookerquartic

function dbookerquartic(ea, M, q, B)
    S, C = sincos(ea)
    C² = C^2

    dS = C
    # dC = -S
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
    upgoing(q)

Calculate the absolute angle of `q` in radians from 315°×π/180 on the complex plane. Smaller
values indicate upgoing waves.
"""
function upgoing(q::Complex)
    # -π/4 is 315°
    a = -π/4 - angle(q)

    # angle() is +/-π
    # `a` can never be > 3π/4, but there is a region on the plane where `a` can be < -π
    a < -π && (a += 2π)
    return abs(a)
end

upgoing(q::Real) = oftype(q, π/4)

"""
    sortquarticroots!(q)

Sort array of quartic roots `q` in place such that the first two correspond to upgoing waves
and the latter two correspond to downgoing waves.

General locations of the four roots of the Booker quartic on the complex plane corresponding
to the:

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

Based on [Pitteway1965] fig. 5.

# References

[Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields, reflexion
    coefficients and polarizations for long radio waves in the lower ionosphere. I.,” Phil.
    Trans. R. Soc. Lond. A, vol. 257, no. 1079, pp. 219–241, Mar. 1965,
    doi: 10.1098/rsta.1965.0004.
"""
function sortquarticroots!(q)
    # It looks like we could just `sort!(q, by=real)`, but the quadrants of each root are
    # not fixed and the positions in the Argand diagram are just approximate.

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

@doc raw"""
    bookerwavefields(ea, M)
    bookerwavefields(T::TMatrix)

Compute the two-column wavefields matrix `e` from the ionosphere with susceptibility tensor
`M` or `T` matrix for the two upgoing wavefields.

The first column of `e` is the evanescent wave and the second is the travelling wave.

```math
e =
    \begin{pmatrix}
         Ex₁  &  Ex₂ \\
        -Ey₁  & -Ey₂ \\
         Hx₁  &  Hx₂ \\
         Hy₁  &  Hy₂
    \end{pmatrix}
```

This function solves the eigenvalue problem ``Te = qe``. First, the Booker quartic is solved
for the roots `q`. Then they are sorted so that the roots associated with the two upgoing
waves can be selected. `e` is solved as the eigenvectors for the two `q`s. An
analytical solution is used where `e[2,:] = 1`.
"""
bookerwavefields

function bookerwavefields(T::TMatrix)
    q, _ = bookerquartic(T)
    sortquarticroots!(q)
    return bookerwavefields(T, q)
end

function bookerwavefields(ea, M)
    q, _ = bookerquartic(ea, M)
    sortquarticroots!(q)
    T = tmatrix(ea, M)
    return bookerwavefields(T, q)
end

function bookerwavefields(T::TMatrix, q)
    # Precompute
    T14T41 = T[1,4]*T[4,1]
    T14T42 = T[1,4]*T[4,2]
    T12T41 = T[1,2]*T[4,1]

    e = MArray{Tuple{4,2},ComplexF64,2,8}(undef)
    @inbounds for i = 1:2
        d = T14T41 - (T[1,1] - q[i])*(T[4,4] - q[i])

        e[1,i] = (T[1,2]*(T[4,4] - q[i]) - T14T42)/d
        e[2,i] = 1
        e[3,i] = q[i]
        e[4,i] = (-T12T41 + T[4,2]*(T[1,1] - q[i]))/d
    end

    # By returning as SArray instead of MArray, the MArray doesn't get hit by GC
    return SArray(e)
end

"""
    bookerwavefields(ea, M, ::Dθ)
    bookerwavefields(T::TMatrix, dT, ::Dθ)

Compute the two-column wavefields matrix `e` as well as its derivative with respect to
``θ``, returned as a tuple `(e, de)` for the ionosphere with susceptibility tensor `M` or
`T` matrix and its derivative with respect to ``θ``, `dT`.
"""
function bookerwavefields(T::TMatrix, dT, ::Dθ)
    q, B = bookerquartic(T)
    sortquarticroots!(q)
    dq = dbookerquartic(T, dT, q, B)
    return bookerwavefields(T, dT, q, dq)
end

function bookerwavefields(ea, M, ::Dθ)
    q, B = bookerquartic(ea, M)
    sortquarticroots!(q)
    dq = dbookerquartic(ea, M, q, B)
    T = tmatrix(ea, M)
    dT = dtmatrix(ea, M)

    return bookerwavefields(T, dT, q, dq)
end

function bookerwavefields(T::TMatrix, dT, q, dq)
    # Precompute
    T14T41 = T[1,4]*T[4,1]
    T14T42 = T[1,4]*T[4,2]
    T12T41 = T[1,2]*T[4,1]

    e = MArray{Tuple{4,2},ComplexF64,2,8}(undef)
    de = similar(e)
    @inbounds for i = 1:2
        den = T14T41 - (T[1,1] - q[i])*(T[4,4] - q[i])
        den² = abs2(den)
        dden = dT[1,4]*T[4,1] - (T[1,1]*dT[4,4] + dT[1,1]*T[4,4]) +
            q[i]*(dT[1,1] + dT[4,4] - 2*dq[i]) + dq[i]*(T[1,1] + T[4,4])

        e1num = T[1,2]*(T[4,4] - q[i]) - T14T42
        e4num = -T12T41 + T[4,2]*(T[1,1] - q[i])
        e[1,i] = e1num/den
        e[2,i] = 1
        e[3,i] = q[i]
        e[4,i] = e4num/den

        # some dT terms == 0
        de1num = T[1,2]*(dT[4,4] - dq[i]) + dT[1,2]*(T[4,4] - q[i]) - dT[1,4]*T[4,2]
        de4num = -dT[1,2]*T[4,1] + T[4,2]*(dT[1,1] - dq[i])
        de[1,i] = (de1num*den - e1num*dden)/den²
        de[2,i] = 0
        de[3,i] = dq[i]
        de[4,i] = (de4num*den - e4num*dden)/den²
    end

    # By returning as SArray instead of MArray, the MArray doesn't get hit by GC
    return SArray(e), SArray(de)
end

@doc raw"""
    bookerreflection(ea, M::SMatrix{3,3})
    bookerreflection(ea, e)

Compute the ionosphere reflection coefficient matrix for a sharply bounded ionosphere from
4×2 wavefields matrix `e` or the susceptibility matrix `M`.

The ionosphere reflection coefficient matrix is computed from a ratio of the downgoing to
upgoing plane waves in the free space beneath the ionosphere [Budden1988] pg. 307.
These are obtained from the two upgoing characteristic waves found from the Booker quartic.
Each make up a column of `e`.

```math
R =
    \begin{pmatrix}
        Ce₁[4] - e₁[1] &  Ce₂[4] - e₂[1] \\
       -Ce₁[2] + e₁[3] & -Ce₂[2] + e₂[3]
    \end{pmatrix}
    \begin{pmatrix}
        Ce₁[4] + e₁[1] &  Ce₂[4] + e₂[1] \\
       -Ce₁[2] - e₁[3] & -Ce₂[2] - e₂[3]
    \end{pmatrix}^{-1}
```

The reflection coefficient matrix for the sharply bounded case is commonly used as a
starting solution for integration of the reflection coefficient matrix through the
ionosphere.

# References

[Budden1988]: K. G. Budden, “The propagation of radio waves: the theory of radio
    waves of low power in the ionosphere and magnetosphere,” First paperback edition.
    New York: Cambridge University Press, 1988.

# Extended help

The set of horizontal field components ``e = (Ex, -Ey, Z₀Hx, Z₀Hy)ᵀ`` can be separated into
an upgoing and downgoing wave, each of which is generally elliptically polarized. A ratio of
the amplitudes of these two waves give a reflection coefficient, except it would only apply
for an incident wave of that particular elliptical polarization. However, the first set of
fields can be linearly combined with a second independent solution for the fields, which
will generally have a different elliptical polarization than the first. Two linear
combinations of the two sets of fields are formed with unit amplitude, linearly polarized
incident waves. The reflected waves then give the components ``R₁₁``, ``R₂₁`` or ``R₁₂``,
``R₂₂`` for the incident wave in the plane of incidence and perpendicular to it,
respectively [Budden1988] pg 552.

The process for determining the reflection coefficient requires resolving the two sets of
fields ``e₁`` and ``e₂`` into the four linearly polarized vacuum modes. The layer of vacuum
can be assumed to be so thin that it does not affect the fields. There will be two upgoing
waves and two downgoing waves, each which has one ``E`` and one ``H`` in the plane of
incidence. If ``f₁, f₂, f₃, f₄`` are the complex amplitudes of the four component waves,
then in matrix notation ``e = Lf`` where ``L`` is the appropriate transformation matrix.

For ``e₁`` and ``e₂``, we can find the corresponding vectors ``f1`` and ``f2`` by
``f1 = L⁻¹e₁``, ``f2 = L⁻¹e₂`` where the two column vectors are partitioned such that
``f1 = (u1, d1)ᵀ`` and ``f2 = (u2, d2)ᵀ`` for upgoing and downgoing 2-element vectors ``u``
and ``d``. From the definition of the reflection coefficient ``R``, ``d = Ru``. Letting
``U = (u1, u2)``, ``D = (d1, d2)``, then ``D = RU`` and the reflection coefficient is
``R = DU⁻¹``. Because the reflection coefficient matrix is a ratio of fields, either ``e₁``
and/or ``e₂`` can be independently multiplied by an arbitrary constant and the value of
``R`` is unaffected.

This function directly computes ``D`` and ``U`` and solves for ``R`` using the right
division operator `R = D/U`.

For additional details, see [Budden1988], chapter 18, section 7.
"""
bookerreflection

function bookerreflection(ea, e)
    C = cos(ea)

    # The order of the two upgoing waves doesn't matter. Swapping the first and second
    # columns of `e` will result in the same reflection coefficient matrix.

    D = SMatrix{2,2}(C*e[4,1]-e[1,1], -C*e[2,1]+e[3,1], C*e[4,2]-e[1,2], -C*e[2,2]+e[3,2])
    U = SMatrix{2,2}(C*e[4,1]+e[1,1], -C*e[2,1]-e[3,1], C*e[4,2]+e[1,2], -C*e[2,2]-e[3,2])

    R = D/U

    return R
end

function bookerreflection(ea, M::SMatrix{3,3})
    e = bookerwavefields(ea, M)
    return bookerreflection(ea, e)
end

"""
    bookerreflection(ea, M, ::Dθ)

Compute the ionosphere reflection coefficient matrix ``R`` for a sharply bounded
ionosphere with susceptibility tensor `M`, as well as its derivative ``dR/dθ`` returned as
the tuple `(R, dR)`.
"""
function bookerreflection(ea, M, ::Dθ)
    S, C = sincos(ea)

    e, de = bookerwavefields(ea, M, Dθ())

    D = SMatrix{2,2}(C*e[4,1]-e[1,1], -C*e[2,1]+e[3,1], C*e[4,2]-e[1,2], -C*e[2,2]+e[3,2])
    dD = SMatrix{2,2}(-S*e[4,1] + C*de[4,1] - de[1,1], S*e[2,1] - C*de[2,1] + de[3,1],
                      -S*e[4,2] + C*de[4,2] - de[1,2], S*e[2,2] - C*de[2,2] + de[3,2])

    U = SMatrix{2,2}(C*e[4,1]+e[1,1], -C*e[2,1]-e[3,1], C*e[4,2]+e[1,2], -C*e[2,2]-e[3,2])
    dU = SMatrix{2,2}(-S*e[4,1] + C*de[4,1] + de[1,1], S*e[2,1] - C*de[2,1] - de[3,1],
                      -S*e[4,2] + C*de[4,2] + de[1,2], S*e[2,2] - C*de[2,2] - de[3,2])

    R = D/U
    dR = dD/U + D*(-U\dU/U)

    return R, dR
end
