"""
    EigenAngle{T}

Waveguide EigenAngle `θ` (rad), as well as automatically computed values:

  - `cosθ`
  - `sinθ`
  - `cos²θ`
  - `sin²θ`

Technically this is an angle of incidence from the vertical, and not necessarily an
*eigen* angle unless it is found to be associated with a propagated mode in the waveguide.
"""
struct EigenAngle{T}
    θ::T  # radians, because df/dθ are in radians
    cosθ::T
    sinθ::T
    cos²θ::T
    sin²θ::T

    function EigenAngle{T}(θ::T) where T <: Number
        rθ, iθ = reim(θ)
        ((abs(rθ) > 2π) || (abs(iθ) > 2π)) && @warn "θ should be in radians"
        S, C = sincos(θ)
        C² = C^2
        S² = 1 - C²
        new(θ, C, S, C², S²)
    end
end
EigenAngle(θ::T) where T <: Number = EigenAngle{T}(θ)
EigenAngle(ea::EigenAngle) = ea

Base.eltype(::Type{EigenAngle{T}}) where {T} = T

function Base.isless(lhs::EigenAngle{T}, rhs::EigenAngle{T}) where {T <: Complex}
    # Compare complex EigenAngle by real and then imag
    # By defining `isless()`, we get the wanted sorting of `EigenAngle`s
    return isless((real(lhs.θ), imag(lhs.θ)), (real(rhs.θ), imag(rhs.θ)))
end

# function Base.show(io::IO, ::MIME"text/plain", e::EigenAngle{T}) where {T}
#     println(io, "EigenAngle{$T}: ")
#     print(io, " ", getfield(e, :θ))
# end

function Base.show(io::IO, e::EigenAngle{T}) where {T}
    compact = get(io, :compact, false)

    if compact
        print(io, e.θ)
    else
        print(io, e.θ)
    end

end

function Base.show(io::IO, ::MIME"text/plain", e::EigenAngle{T}) where {T}
    println(io, "EigenAngle{$T}: ")
    show(io, e)
end
