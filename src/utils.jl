#==
Utility functions
==#

"""
    amplitude(e)

Compute field amplitude in dB.
"""
amplitude(e) = 10log10(abs2(e))  # == 20log10(abs(E))

"""
    amplitudephase(e)

Compute field amplitude in dB and phase in radians and return as (`amplitude`, `phase`).
"""
function amplitudephase(e::Number)
    a = amplitude(e)
    p = angle(e)

    return a, p
end

function amplitudephase(e)
    a = similar(e, Float64)
    p = similar(a)
    for (c, ec) in pairs(eachcol(e))
        for i in eachindex(ec)
            a[i,c] = amplitude(ec[i])
            p[i,c] = angle(ec[i])
        end
    end
    return a, p
end

"""
    unwrap!(x)

Unwrap a phase vector `x` in radians in-place.
"""
unwrap!

function unwrap!(x)
	for ec in eachcol(x)
        v = first(ec)
        for k in eachindex(ec)
            if isfinite(v)
                ec[k] = v = v + rem2pi(ec[k]-v, RoundNearest)
            end
        end
    end
	return x
end

unwrap!(x::Number) = x

"""
    pow23(x)

Efficiently compute ``x^(2/3)``.
"""
pow23

pow23(x::Real) = cbrt(x)^2
pow23(z::Complex) = exp(2/3*log(z))
