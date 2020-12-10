#==
Functions related to the Booker quartic equation and the four waves associated with the
quartic roots.

Besides calculating the quartic roots themselves, this file also includes functions
related to calculating wavefields and vacuum reflection coefficients from the quartic roots.
==#

"""
    bookerquartic

Compute roots `q` and the coefficients `B` of the Booker quartic.

See also: [`dbookerquartic`](@ref)
"""
function bookerquartic end

"`bookerquartic(ea::EigenAngle, M)`"
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

"`bookerquartic(T::TMatrix)`"
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
    dbookerquartic

Compute derivative `dq` of the Booker quartic roots with respect to ``θ``.

See also: [`bookerquartic`](@ref)
"""
function dbookerquartic end

"`dbookerquartic(ea::EigenAngle, M, q, B)`"
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

"`dbookerquartic(T::TMatrix, dT, q, B)`"
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

Based on [^Pitteway1965] fig. 5.

[^Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields, reflexion
    coefficients and polarizations for long radio waves in the lower ionosphere. I.,” Phil.
    Trans. R. Soc. Lond. A, vol. 257, no. 1079, pp. 219–241, Mar. 1965,
    doi: 10.1098/rsta.1965.0004.
"""
function sortquarticroots!(q)
    # It looks like we could just `sort!(q, by=real)`, but the quadrants of each root are
    # not fixed and the positions in the Argand diagram are just approximate.

    length(q) == 4 || @warn "length of `q` is not 4"

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
    bookerwavefields

Compute the two-column wavefields matrix `e` from the Booker quartic for the two upgoing
wavefields where subscript `1` is the evanescent wave and `2` is the travelling wave.

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
function bookerwavefields end

"`bookerwavefields(T::TMatrix)`"
function bookerwavefields(T::TMatrix)
    q, B = bookerquartic(T)
    sortquarticroots!(q)
    return bookerwavefields(T, q)
end

"`bookerwavefields(ea::EigenAngle, M)`"
function bookerwavefields(ea::EigenAngle, M)
    q, B = bookerquartic(ea, M)
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
    bookerwavefields(T::TMatrix, dT, ::Dθ)

Compute the two-column wavefields `e` as well as its derivative with respect to ``θ``,
returned as a tuple `e, de`.
"""
function bookerwavefields(T::TMatrix, dT, ::Dθ)
    q, B = bookerquartic(T)
    sortquarticroots!(q)
    dq = dbookerquartic(T, dT, q, B)
    return bookerwavefields(T, dT, q, dq)
end

"`bookerwavefields(ea::EigenAngle, M, ::Dθ)`"
function bookerwavefields(ea::EigenAngle, M, ::Dθ)
    q, B = bookerquartic(ea, M)
    sortquarticroots!(q)
    dq = dbookerquartic(ea, M, q, B)
    T = tmatrix(ea, M)
    dT = dtmatrix(ea, M)

    return bookerwavefields(T, dT, q, dq)
end

"`bookerwavefields(T::TMatrix, dT, q, dq)`"
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
    bookerreflection

The reflection coefficient matrix is calculated from a ratio of the downgoing to upgoing
plane waves in the free space beneath the ionosphere [^Budden1988] pg. 307. These are
obtained from the two upgoing characteristic waves found from the Booker quartic. Each make
up a column of `e`.

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

See also: [`bookerwavefields`](@ref)

# References

[^Budden1988]: K. G. Budden, “The propagation of radio waves,”
    Cambridge University Press, 1988.
"""
function bookerreflection end

"`bookerreflection(ea::EigenAngle, e)`"
function bookerreflection(ea::EigenAngle, e)
    C = ea.cosθ

    # The order of the two upgoing waves doesn't matter. Swapping the first and second
    # columns of `e` will result in the same reflection coefficient matrix.

    D = SMatrix{2,2}(C*e[4,1]-e[1,1], -C*e[2,1]+e[3,1], C*e[4,2]-e[1,2], -C*e[2,2]+e[3,2])
    U = SMatrix{2,2}(C*e[4,1]+e[1,1], -C*e[2,1]-e[3,1], C*e[4,2]+e[1,2], -C*e[2,2]-e[3,2])

    R = D/U

    return R
end

"`bookerreflection(ea::EigenAngle, M::SMatrix{3,3})`"
function bookerreflection(ea::EigenAngle, M::SMatrix{3,3})
    e = bookerwavefields(ea, M)
    return bookerreflection(ea, e)
end

"""
    bookerreflection(ea::EigenAngle, M, ::Dθ)

Compute the ionosphere reflection coefficient matrix ``R`` for a sharply bounded
ionosphere as well as its derivative ``dR/dθ`` returned as tuple `R, dR`.
"""
function bookerreflection(ea::EigenAngle, M, ::Dθ)
    S, C = ea.sinθ, ea.cosθ

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
