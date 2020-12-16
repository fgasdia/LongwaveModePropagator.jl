#==
Utility functions
==#

"""
    unwrap!(x)

Unwrap a phase vector `x` in radians in-place.
"""
function unwrap!(x)
	v = first(x)
	@inbounds for k in eachindex(x)
		x[k] = v = v + rem2pi(x[k]-v, RoundNearest)
	end
	return x
end

"""
    pow23(x)

Efficiently compute ``x^(2/3)``.
"""
function pow23 end
pow23(x::Real) = cbrt(x)^2
pow23(z::Complex) = exp(2/3*log(z))
