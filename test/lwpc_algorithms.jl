
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
