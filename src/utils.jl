#==
Utility functions
==#

"""
    unwrap!(x)

Unwrap a phase vector `x` in radians in-place.
"""
unwrap!

function unwrap!(x::AbstractVector)
	v = first(x)
	@inbounds for k in eachindex(x)
		x[k] = v = v + rem2pi(x[k]-v, RoundNearest)
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
