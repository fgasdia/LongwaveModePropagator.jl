using StaticArrays
using PolynomialRoots

"""
Solution for the roots of a fourth-order polynomial (quartic equation) taken from Burnside
and Panton (1904), the Theory of Equations. A summary of the equations and process is given
in Radio Science Vol. 3, Aug 1968, pg. 792-795.
"""
function solvequartic(B4, B3, B2, B1, B0)
    # More suitable form for solution:
    # q⁴ + 4b₃q³ + 6b₂q² + 4b₁q + b₀ = 0
    b3 = B3/4B4
    b2 = B2/6B4
    b1 = B1/4B4
    b0 = B0/B4

    H = b2 - b3^2
    I = b0 - 4*b3*b1 + 3b2^2
    G = b1 - 3*b3*b2 + 2b3^3
    h = -I/12
    g = -G^2/4 - H*(H^2 + 3h)
    p = (-g + sqrt(g^2 + 4h^3))/2
    magpos = abs(real(p)) + abs(imag(p))
    ppos = p
    p = (-g - sqrt(g^2 + 4h^3))/2
    magneg = abs(real(p)) + abs(imag(p))
    magpos > magneg && (p = ppos)

    ω2 = complex(-0.5, sqrt(3)/2)
    ω3 = complex(-0.5, -sqrt(3)/2)
    s1 = exp(log(p)/3)  # cube root of p
    s2 = ω2*s1
    s3 = ω3*s1

    r1 = sqrt(s1 - h/s1 - H)
    r2 = sqrt(s2 - h/s2 - H)
    r3 = sqrt(s3 - h/s3 - H)

    if abs(G) >= 1e-20
        real(-2r1*r2*r3/G) ≈ -1 && (r3 = -r3)
    end

    q = @MVector [r1+r2+r3, r1-r2-r3, -r1+r2-r3, -r1-r2+r3]
    q .-= b3

    # First order Newton-Raphson iterative improvement
    for j = 1:4
        dlqmin = 9.9e9
        for jj = 1:4
            if jj != j
                dlq = abs(q[j] - q[jj])
                dlq < dlqmin && (dqlmin = dlq)
            end
        end
        dlqmax = dlqmin/3

        lastiter = false
        ncount = 1
        delq = 0.0
        while !lastiter & (ncount <= 10)
            f = (((B4*q[j] + B3)*q[j] + B2)*q[j] + B1)*q[j] + B0
            dfdq = ((4B4*q[j] + 3B3)*q[j] + 2B2)*q[j] + B1
            delq = -f/dfdq
            abs(delq) > dlqmax && (delq *= delqmax/abs(delq))
            q[j] += delq
            abs(delq/q[j]) < 1e-4 && (lastiter = true)
            ncount += 1
        end
        if abs(delq/q[j]) > 1e-2
            error("q fails to converge")
        end
    end
    return q
end


"""
Solution for the roots of a fourth-order polynomial (quartic equation) taken from Burnside
and Panton (1904), the Theory of Equations. A summary of the equations and process is given
in Radio Science Vol. 3, Aug 1968, pg. 792-795.
"""
function solvequartic(B4, B3, B2, B1, B0)
    # More suitable form for solution:
    # q⁴ + 4b₃q³ + 6b₂q² + 4b₁q + b₀ = 0
    b3 = B3/4B4
    b2 = B2/6B4
    b1 = B1/4B4
    b0 = B0/B4

    H = b2 - b3^2
    I = b0 - 4*b3*b1 + 3b2^2
    G = b1 - 3*b3*b2 + 2b3^3
    h = -I/12
    g = -G^2/4 - H*(H^2 + 3h)
    p = (-g + sqrt(g^2 + 4h^3))/2
    magpos = abs(real(p)) + abs(imag(p))
    ppos = p
    p = (-g - sqrt(g^2 + 4h^3))/2
    magneg = abs(real(p)) + abs(imag(p))
    magpos > magneg && (p = ppos)

    ω2 = complex(-0.5, sqrt(3)/2)
    ω3 = complex(-0.5, -sqrt(3)/2)
    s1 = exp(log(p)/3)  # cube root of p
    s2 = ω2*s1
    s3 = ω3*s1

    r1 = sqrt(s1 - h/s1 - H)
    r2 = sqrt(s2 - h/s2 - H)
    r3 = sqrt(s3 - h/s3 - H)

    if abs(G) >= 1e-20
        real(-2r1*r2*r3/G) ≈ -1 && (r3 = -r3)
    end

    q = @MVector [r1+r2+r3, r1-r2-r3, -r1+r2-r3, -r1-r2+r3]
    q .-= b3

    # First order Newton-Raphson iterative improvement
    for j = 1:4
        dlqmin = 9.9e9
        for jj = 1:4
            if jj != j
                dlq = abs(q[j] - q[jj])
                dlq < dlqmin && (dqlmin = dlq)
            end
        end
        dlqmax = dlqmin/3

        lastiter = false
        ncount = 1
        delq = 0.0
        while !lastiter & (ncount <= 10)
            f = (((B4*q[j] + B3)*q[j] + B2)*q[j] + B1)*q[j] + B0
            dfdq = ((4B4*q[j] + 3B3)*q[j] + 2B2)*q[j] + B1
            delq = -f/dfdq
            abs(delq) > dlqmax && (delq *= delqmax/abs(delq))
            q[j] += delq
            abs(delq/q[j]) < 1e-4 && (lastiter = true)
            ncount += 1
        end
        if abs(delq/q[j]) > 1e-2
            error("q fails to converge")
        end
    end
    return q
end

test1 = [complex(1.31, 0.21), complex(0.42, 0.13), complex(32.1, 3.1), complex(0.13, 0.4), complex(2.1, 0.9)]
test2 = [complex(2.31, 1.21), complex(0.92, 0.3), complex(23.1, 4.1), complex(0.35, 0.4), complex(1.1, 1.9)]

test11 = [1, test1[2]/(4test1[1]), test1[3]/(6test1[1]), test1[4]/(4test1[1]), test1[5]/test1[1]]
# [4, 3, 2, 1, 0]
sol1 = solvequartic(test1[1], test1[2], test1[3], test1[4], test1[5])
# sol1 = solvequartic(test11[1], test11[2], test11[3], test11[4], test11[5])

# [0, 1, 2, 3, 4]
rtest1 = reverse(test1)
sol2 = roots(rtest1)  # gives same results in less time/memory!

srtest1 = SVector(test1[5], test1[4], test1[3], test1[2], test1[1])
sol3 = roots(srtest1)


q = copy(sol2)

function jnk(q)
l = 0
for m = 2:4, n = m:4
    if imag(q[n]) <= 0
        l += 1
        q[n], q[m-1] = q[m-1], q[n]
    end
end
return l
end

sort!(q, by=imag)


anglediff = @MVector [0., 0., 0., 0.]
if l != 2
    for n = 1:4
        angq = rad2deg(angle(q[n]))
        angq < 0 && (angq += 360)
        angq < 135 && (angq += 360)
        anglediff[n] = abs(angq - 315)
    end
    for nm = 2:4, n = nm:4
        if anglediff[n] <= anglediff[nm-1]
            anglediff[n], anglediff[nm-1] = anglediff[nm-1], anglediff[n]
            q[n], q[nm-1] = q[nm-1], q[n]
        end
    end
end

function angdiff(qval)
    angq = rad2deg(angle(qval))
    angq < 0 && (angq += 360)
    angq < 135 && (angq += 360)
    abs(angq - 315)
end




# Test if R equations are the same (without Δ since they both use it in the same way)
T1 = complex(0.987, 0.013)
T2 = complex(0.565, 0.12)
P1 = complex(1.31, 0.7)
P2 = complex(-1.6, 0.003)
q1 = complex(1.22, -0.019)
q2 = complex(0.74, -0.108)
C = cos(complex(0.672, -0.12))

Rpar = ((T1*C - P1)*(C + q2) - (T2*C - P2)*(C + q1))
Xpar = (T1*(C+q2) - T2*(C+q1))*2




"""
Return a cube roots of a complex number `z`.

See https://math.stackexchange.com/questions/394432/cube-roots-of-the-complex-numbers-1i
"""
function complexcbrts(z::Complex)
    r = abs(z)
    rroot = cbrt(r)
    θ = angle(z)
    dθ = θ/3

    r1 = rroot*cis(dθ)
    r2 = rroot*cis(2π/3 + dθ)
    r3 = rroot*cis(2π*2/3 + dθ)

    return r1, r2, r3
end

function bookerquartic(B4, B3, B2, B1, B0)
    b3 = B3/(4*B4)
    b2 = B2/(6*B4)
    b1 = B1/(4*B4)
    b0 = B0/B4

    b3² = b3^2

    H = b2 - b3²
    Iv = b0 - 4*b3*b1 + 3*b2^2
    G = b1 - 3*b3*b2 + 2*b3²
    h = -Iv/12
    g = -G^2/4 - H*(H^2 + 3*h)

    tmpsqrt = sqrt(complex(g^2 + 4*h^3))
    p = (-g + tmpsqrt)/2
    abs(p) < 1e-10 && (p = (-g - tmpsqrt)/2)

    # LWPC hard codes two of the unit cube roots and uses exp(log(p)/3) for cbrt
    s1, s2, s3 = complexcbrts(p)

    r1 = sqrt(s1 - h/s1 - H)
    r2 = sqrt(s2 - h/s2 - H)
    r3 = sqrt(s3 - h/s3 - H)

    # If G != 0, G should = +1, if G = -1, switch sign of one r
    abs2(G) > 0 && real(-2*r1*r2*r3/G) < 0 && (r1 = -r1)

    q1 = r1 + r2 + r3 - b3
    q2 = r1 - r2 - r3 - b3
    q3 = -r1 + r2 - r3 - b3
    q4 = -r1 - r2 + r3 - b3

    # LWPC then does first order Newton-Raphson iterative improvement, probably
    # because the above algorithm is AWFUL in comparison to the above algorithms
    # It's not clear if I need to an iterative refinement to the above closed
    # form solution

    return q1, q2, q3, q4
end
