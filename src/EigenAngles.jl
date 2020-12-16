"""
    EigenAngle

Plane-wave propagation direction in angle `θ` from the vertical in radians.

Common trigonometric function values of this angle are calculated to increase performance.

# Fields

- `θ::ComplexF64`
- `cosθ::ComplexF64`
- `sinθ::ComplexF64`
- `secθ::ComplexF64`
- `cos²θ::ComplexF64`
- `sin²θ::ComplexF64`

`Real` `θ` will be automatically converted to `Complex`.

!!! note

    Technically this is an angle of incidence from the vertical, and not necessarily an
    _eigen_angle unless it is found to be associated with a propagated mode in the waveguide.
"""
struct EigenAngle
    θ::ComplexF64  # radians, because df/dθ are in radians
    cosθ::ComplexF64
    sinθ::ComplexF64
    secθ::ComplexF64  # == 1/cosθ
    cos²θ::ComplexF64
    sin²θ::ComplexF64
end

EigenAngle(ea::EigenAngle) = ea
EigenAngle(θ::Real) = EigenAngle(complex(θ))
function EigenAngle(θ::Complex)
    rθ, iθ = reim(θ)
    ((abs(rθ) > 2π) || (abs(iθ) > 2π)) && @warn "θ > 2π. Make sure θ is in radians."
    S, C = sincos(θ)
    Cinv = inv(C)  # == sec(θ)
    C² = C^2
    S² = 1 - C²
    EigenAngle(θ, C, S, Cinv, C², S²)
end

"""
    isless(x::EigenAngle, y::EigenAngle)

Calls `isless` for complex `EigenAngle` by `real` and then `imag` component.

By defining `isless`, calling `sort` on `EigenAngle`s sorts by real component first and
imaginary component second.
"""
function Base.isless(x::EigenAngle, y::EigenAngle)
    return isless((real(x.θ), imag(x.θ)), (real(y.θ), imag(y.θ)))
end

function Base.show(io::IO, e::EigenAngle)
    compact = get(io, :compact, false)

    if compact
        print(io, e.θ)
    else
        print(io, e.θ)
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", e::EigenAngle)
    println(io, "EigenAngle: ")
    show(io, e)
end

"""
    isdetached(ea::EigenAngle, frequency::Frequency; params=LMPParams())

Return `true` if this is likely an earth detached (whispering gallery) mode according to the
criteria in [^Pappert1981] eq. 1 with the additional criteria that the frequency be above
100 kHz.

# References

[^Pappert1981]: R. A. Pappert, “LF daytime earth ionosphere waveguide calculations,” Naval
    Ocean Systems Center, San Diego, CA, NOSC/TR-647, Jan. 1981.
    [Online]. Available: https://apps.dtic.mil/docs/citations/ADA096098.
"""
function isdetached(ea::EigenAngle, frequency::Frequency; params=LMPParams())
    C, C² = ea.cosθ, ea.cos²θ
    f, k = frequency.f, frequency.k
    @unpack earthradius, curvatureheight = params
    α = 2/earthradius

    return f > 100e3 && real(2im/3*(k/α)*(C²- α*curvatureheight)^(3/2)) > 12.4
end

"""
    referencetoground(ea::EigenAngle; params=LMPParams())

Reference `ea` from `params.curvatureheight` to ground (``z = 0``).
"""
function referencetoground(ea::EigenAngle; params=LMPParams())
    # see, e.g. PS71 pg 11
    @unpack earthradius, curvatureheight = params
    return EigenAngle(asin(ea.sinθ/sqrt(1 - 2/earthradius*curvatureheight)))
end

"""
    referencetoground(sinθ; params=LMPParams())

Reference ``sin(θ)`` from `params.curvatureheight` to the ground and return ``sin(θ₀)`` at
the ground.
"""
function referencetoground(sinθ; params=LMPParams())
    @unpack earthradius, curvatureheight = params
    return sinθ/sqrt(1 - 2/earthradius*curvatureheight)
end

"""
    attenuation(ea, frequency::Frequency; params=LMPParams())

Compute attenuation of `ea` at the ground.

This function internally references `ea` to the ground.
"""
function attenuation(ea, frequency::Frequency; params=LMPParams())
    ea = EigenAngle(ea)
    S₀ = referencetoground(ea.sinθ)
    neper2dB = 20log10(exp(1))  # 1 Np ≈ 8.685 dB
    return -neper2dB*frequency.k*imag(S₀)*1e6
end

"""
    phasevelocity(ea)

Compute the relative phase velocity ``v/c``.

This function internally references `ea` to the ground.
"""
function phasevelocity(ea)
    ea = EigenAngle(ea)
    S₀ = referencetoground(ea.sinθ)
    return 1/real(S₀)
end
