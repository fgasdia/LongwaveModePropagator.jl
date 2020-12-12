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
    referencetoground(ea::EigenAngle; params=LMPParams())

Reference `ea` from `params.curvatureheight` to ground (``z = 0``).
"""
function referencetoground(ea::EigenAngle; params=LMPParams())
    # see, e.g. PS71 pg 11
    @unpack earthradius, curvatureheight = params
    return EigenAngle(asin(ea.sinθ/sqrt(1 - 2/earthradius*curvatureheight)))
end

"""
    referencetoground(x; params=LMPParams())

Reference ``x == sin(ea)`` from `params.curvatureheight` to the ground and return `sin(ea)` at
the ground.
"""
function referencetoground(x; params=LMPParams())
    @unpack earthradius, curvatureheight = params
    return x/sqrt(1 - 2/earthradius*curvatureheight)
end
