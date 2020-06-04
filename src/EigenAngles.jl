"""
    EigenAngle

Waveguide EigenAngle `θ` (rad), as well as automatically computed values:

  - `cosθ`
  - `sinθ`
  - `cos²θ`
  - `sin²θ`

Technically this is an angle of incidence from the vertical, and not necessarily an
*eigen* angle unless it is found to be associated with a propagated mode in the waveguide.
"""
struct EigenAngle
    θ::ComplexF64  # radians, because df/dθ are in radians
    cosθ::ComplexF64
    sinθ::ComplexF64
    secθ::ComplexF64  # == 1/cosθ
    cos²θ::ComplexF64
    sin²θ::ComplexF64

    function EigenAngle(θ::ComplexF64)
        rθ, iθ = reim(θ)
        ((abs(rθ) > 2π) || (abs(iθ) > 2π)) && @warn "θ > 2π. Make sure θ is in radians."
        S, C = sincos(θ)
        Cinv = inv(C)  # == sec(θ)
        C² = C^2
        S² = 1 - C²
        new(θ, C, S, Cinv, C², S²)
    end
end
EigenAngle(θ::Real) = EigenAngle(complex(float(θ)))
EigenAngle(ea::EigenAngle) = ea
Base.eltype(::Type{EigenAngle}) = ComplexF64

function Base.isless(lhs::EigenAngle, rhs::EigenAngle)
    # Compare complex EigenAngle by real and then imag
    # By defining `isless()`, we get the wanted sorting of `EigenAngle`s
    return isless((real(lhs.θ), imag(lhs.θ)), (real(rhs.θ), imag(rhs.θ)))
end

# function Base.show(io::IO, ::MIME"text/plain", e::EigenAngle{T}) where {T}
#     println(io, "EigenAngle{$T}: ")
#     print(io, " ", getfield(e, :θ))
# end

function Base.show(io::IO, e::EigenAngle)
    compact = get(io, :compact, false)

    if compact
        print(io, e.θ)
    else
        print(io, e.θ)
    end

end

function Base.show(io::IO, ::MIME"text/plain", e::EigenAngle)
    println(io, "EigenAngle{$T}: ")
    show(io, e)
end
