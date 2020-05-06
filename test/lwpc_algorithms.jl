
function what(v)
    angq = angle(v)
    imag(v) < 0 && (angq += 2π)
    angq < -π/4 && (angq += 2π)
    return angq
end

function what2(v)
    angq = rad2deg(angle(v))
    imag(v) < 0 && (angq += 360)
    angq < 315 && (angq += 360)
    return angq
end

function lwsort(v)
    angq = rad2deg(angle(v))
    angq < 0 && (angq += 360)
    angq < 135 && (angq += 360)
    return abs(angq-315)
end


function lwsortfull(v)
    q = copy(v)
    diffq = similar(q, Float64)

    l = 0
    for m = 2:4
        for n = m:4
            if imag(q[n]) <= 0
                l += 1
                q[n], q[m-1] = q[m-1], q[n]
            end
        end
    end
    if l != 2
        for n = 1:4
            diffq[n] = lwsort(q[n])
        end
        for nm = 2:4
            for n = nm:4
                if diffq[n] <= diffq[nm-1]
                    diffq[n], diffq[nm-1] = diffq[nm-1], diffq[n]
                    q[n], q[nm-1] = q[nm-1], q[n]
                end
            end
        end
    end
    return q
end

function wfsortfull(v)
    q = copy(v)
    diffq = similar(q, Float64)

    imag(q[1]) <= 0 ? l = 1 : l = 0

    m = l + 1
    while m < 4
        n = m + 1
        while n <= 4
            if imag(q[n]) <= 0
                l += 1
                q[n], q[m] = q[m], q[n]
                m += 1
            end
            n+= 1
        end
        m += 1
    end

    if l != 2
        m = 1
        while m <= 4
            diffq[m] = what2(q[m])
            m += 1
        end
        m = 2
        while m <= 4
            n = m
            while n <= 4
                if diffq[n] <= diffq[m-1]
                    diffq[n], diffq[m-1] = diffq[m-1], diffq[n]
                    q[n], q[m-1] = q[m-1], q[n]
                end
                n += 1
            end
            m += 1
        end
    end
    return q
end


function lwpcreflectioncoeffs(ea::EigenAngle, e1, e2)
    # From wf_r_mtrx.for

    C = ea.cosθ

    g12 = e1[1]*e2[2] - e2[1]*e1[2]
    g13 = e1[1]*e2[3] - e2[1]*e1[3]
    g14 = e1[1]*e2[4] - e2[1]*e1[4]
    g23 = e1[2]*e2[3] - e2[2]*e1[3]
    g24 = e1[2]*e2[4] - e2[2]*e1[4]
    g34 = e1[3]*e2[4] - e2[3]*e1[4]

    den = -g13 + C*(g34 - g12 + C*g24)

    d11 = g13 + C*(g34 + g12 + C*g24)
    d22 = g13 + C*(-g34 - g12 + C*g24)
    d12 = 2*C*g14
    d21 = 2*C*g23

    return SMatrix{2,2,eltype(den),4}(d11/den, d21/den, d12/den, d22/den)
end

vacuumR = vacuumreflectioncoeffs(ea, e1[end], e2[end])
@test vacuumR ≈ lwpcreflectioncoeffs(ea, e1[end], e2[end])

function lwpcscale(p1, p2)
    e1, e2 = MVector(p1), MVector(p2)

    # aterm → dot(e1, e1)
    aterm = 0
    for i = 1:4
        aterm += abs2(e1[i])
    end

    # term → dot(e1, e2)
    term = 0
    for i = 1:4
        term += conj(e1[i])*e2[i]
    end

    # term → dot(e1, e2)/dot(e1, e1)
    term /= aterm

    # e2 → e2 - dot(e1, e2)/dot(e1, e1)
    for i = 1:4
        e2[i] -= term*e1[i]
    end

    # bterm → dot(e2, e2)
    bterm = 0
    for i = 1:4
        bterm += abs2(e2[i])
    end

    # Normalize both vectors
    aterm = 1/sqrt(aterm)
    bterm = 1/sqrt(bterm)
    for i = 1:4
        e1[i] *= aterm
        e2[i] *= bterm
    end

    return SVector(e1), SVector(e2)
end
