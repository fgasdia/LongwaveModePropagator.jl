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
    @inbounds for i in eachindex(e)
        a[i] = amplitude(e[i])
        p[i] = angle(e[i])
    end
    return a, p
end

"""
    unwrap!(x)

Unwrap a phase vector `x` in radians in-place.
"""
unwrap!

function unwrap!(x)
	v = first(x)
	@inbounds for k in eachindex(x)
        if isfinite(v)
            x[k] = v = v + rem2pi(x[k]-v, RoundNearest)
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
