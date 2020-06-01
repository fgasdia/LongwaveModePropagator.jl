# CODATA 2018 NIST SP 961
const c₀ = 299792458.0  # speed of light in vacuum, m/s
const μ₀ = 1.25663706212e-6  # vacuum permeability, H/m
const ϵ₀ = 1/(μ₀*c₀^2)  # vacuum permittivity, F/m
const Z₀ = 376.730313668  # vacuum impedance, Ω

const Rₑ = 6366e3  # earth radius, m
const CURVATURE_HEIGHT = 50e3  # reference height for Earth curvature, m lwpm.for defines as 50 km

########

"""
    Species{F, G}

Ionosphere constituent `Species` with a `charge` (C), `mass` (kg),
`numberdensity(z)` (m⁻³), and `collisionfrequency(z)` (s⁻¹).
"""
struct Species{F<:Function, G<:Function}
    charge::Float64  # C
    mass::Float64  # kg
    numberdensity::F  # function that obtains number density (m⁻³)
    collisionfrequency::G  # function that obtains collision frequency (s⁻¹)
end
# TODO: Suggest to user that F (and maybe G) be an interpolated object, which
# may be faster than the calls to `exp` that will likely appear in F. The calls
# to `exp` accounts for ~30% of `susceptibility` runtime. `G` also has a call to
# exp typically, but this is much less time.

########

"""
    BField

Background magnetic field vector (T).

`B` should be in Teslas, angles in radians!
"""
struct BField
    B::Float64
    dcl::Float64
    dcm::Float64
    dcn::Float64
end

"""
    BField(B, dip, azimuth)

Return a `BField` vector of strength `B` (T), `dip` angle (rad) from the horizontal, and
`azimuth` angle (rad) from ???.
"""
function BField(B, dip, azimuth)
    # BUG: How can we avoid this?
    # Look at the Booker quartic roots for dip angles -1,+1. They are way
    # outside the normal. The roots aren't accurately found, although I'm not
    # sure where the issue is.
    abs(rad2deg(dip)) <= 1 && @warn "magnetic dip angles between ±1° have known numerical issues."
    abs(dip) > 2π && @warn "magnetic dip angle should be in radians"
    abs(azimuth) > 2π && @warn "magnetic azimuth angle should be in radians"

    Sdip, Cdip = sincos(dip)
    Saz, Caz = sincos(azimuth)
    BField(B, Cdip*Caz, Cdip*Saz, -Sdip)
end

"""
    dip(b::BField)

Return the dip angle (rad) from a `BField` vector `b`.
"""
dip(b::BField) = -asin(b.dcn)

"""
    azimuth(b::BField)

Return the azimuth angle (rad) from a `BField` vector `b`.
"""
azimuth(b::BField) = atan(b.dcm,b.dcl)


########
# Useful plasma functions

"""
    plasmafrequency(n, q, m)

Return the plasma frequency ``ωₚ = √(nq²/(ϵ₀m))`` for a ``cold'' plasma.
"""
plasmafrequency(n, q, m) = sqrt(n*q^2/(ϵ₀*m))
plasmafrequency(z, s::Species) = plasmafrequency(s.numberdensity(z), s.charge, s.mass)

"""
    gyrofrequency(q, m, B)

Return the _signed_ gyrofrequency ``Ω = qB/m``.
"""
gyrofrequency(q, m, B) = q*B/m
gyrofrequency(s::Species, b::BField) = gyrofrequency(s.charge, s.mass, b.B)

########

"""
    Ground(ϵᵣ, σ)

Isotropic earth ground characterized by a relative permittivity `ϵᵣ` and a
conductivity `σ`.
"""
struct Ground
    ϵᵣ::Int
    σ::Float64
end

########

# NOTE:
# @kristoffer.carlson suggests this is faster than e.g. dispatch
# https://discourse.julialang.org/t/dispatch-and-symbols/21162/5
# However, there are limitations to enum in Julia 1.0, mainly the enum fields are in scope
# everywhere.
# TODO: Use SuperEnum or something else with improved scoping? >
# https://discourse.julialang.org/t/encapsulating-enum-access-via-dot-syntax/11785/9
@enum FieldComponent begin
    FC_Ez
    FC_Ey
    FC_Ex
    FC_ALL
end

########

# TODO: Move EigenAngle out of Geophysics

"""
    EigenAngle{T}

Waveguide EigenAngle `θ` (rad), as well as automatically computed values:

  - `cosθ`
  - `sinθ`
  - `cos²θ`
  - `sin²θ`

Technically this is an angle of incidence from the vertical, and not necessarily an
*eigen* angle unless it is found to be associated with s propagated mode in the waveguide.
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


########

"""
    waitprofile(z, h′, β)

Return the electron number density in m⁻³ at altitude `z` (m) using Wait's exponential profile
with parameters `h′` (km) and `β` (km⁻¹).

``Nₑ = 1.43 × 10¹³ \\exp(-0.15 h') ⋅ \\exp[(β - 0.15)(z/1000 - h')]``

!!! note
    This profile has no altitude cutoffs (it goes to the ground).

See also: [`electroncollisionfrequency`](@ref), [`ioncollisionfrequency`](@ref)

# References
[^Wait1964]: J. R. Wait and K. P. Spies, “Characteristics of the earth-ionosphere waveguide for VLF radio waves,” U.S. National Bureau of Standards, Boulder, CO, Technical Note 300, Dec. 1964.

[^Thomson1993]: N. R. Thomson, “Experimental daytime VLF ionospheric parameters,” Journal of Atmospheric and Terrestrial Physics, vol. 55, no. 2, pp. 173–184, Feb. 1993, doi: 10.1016/0021-9169(93)90122-F.

[^Cummer1998]: S. A. Cummer, U. S. Inan, and T. F. Bell, “Ionospheric D region remote sensing using VLF radio atmospherics,” Radio Science, vol. 33, no. 6, pp. 1781–1792, Nov. 1998, doi: 10.1029/98RS02381.
"""
function waitprofile(z, h′, β)
    # Using form with single `exp` for speed
    return 1.43e13*exp(-0.15h′ - (β-0.15)*(h′ - z/1000))
end

"""
    waitprofile(z, h′, β, cutoff_low)

When `z` is below `cutoff_low` (m), return 0.
"""
function waitprofile(z, h′, β, cutoff_low)
    if z > cutoff_low
        return waitprofile(z, h′, β)
    else
        # waitprofile(z, h′, β) return type is essentially promote_type(Float64, typeof(z))
        return zero(promote_type(Float64, typeof(z)))
    end
end

"""
    electroncollisionfrequency(z)

Return the electron-neutral collision frequency (s⁻¹) at altitude `z` (m) based on Wait's
conductivity profile.

``νₑ(z) = 1.816 × 10¹¹ \\exp(-0.15e-3 z)``

See also: [`waitprofile`](@ref), [`ioncollisionfrequency`](@ref)

# References
[^Wait1964]: J. R. Wait and K. P. Spies, “Characteristics of the earth-ionosphere waveguide for VLF radio waves,” U.S. National Bureau of Standards, Boulder, CO, Technical Note 300, Dec. 1964.

[^Thomson1993]: N. R. Thomson, “Experimental daytime VLF ionospheric parameters,” Journal of Atmospheric and Terrestrial Physics, vol. 55, no. 2, pp. 173–184, Feb. 1993, doi: 10.1016/0021-9169(93)90122-F.

[^Cummer1998]: S. A. Cummer, U. S. Inan, and T. F. Bell, “Ionospheric D region remote sensing using VLF radio atmospherics,” Radio Science, vol. 33, no. 6, pp. 1781–1792, Nov. 1998, doi: 10.1029/98RS02381.
"""
function electroncollisionfrequency(z)
    return 1.816e11*exp(-0.15e-3*z)  # e-3 converts `z` to km
end

"""
    electroncollisionfrequency(z, cutoff_low)

When `z` is below `cutoff_low` (m), return 0.
"""
function electroncollisionfrequency(z, cutoff_low)
    if z > cutoff_low
        return electroncollisionfrequency(z)
    else
        return zero(promote_type(Float64, typeof(z)))
    end
end

"""
    ioncollisionfrequency(z)

Return ion-neutral collision frequency (s⁻¹) at height `z` (m) from [^Morfitt1976].

``νᵢ(z) = 4.54 × 10⁹ \\exp(-0.15e-3 z)``

See also: [`waitprofile`](@ref), [`electroncollisionfrequency`](@ref)

# References
[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.

[^Cummer1998]: S. A. Cummer, U. S. Inan, and T. F. Bell, “Ionospheric D region remote sensing using VLF radio atmospherics,” Radio Science, vol. 33, no. 6, pp. 1781–1792, Nov. 1998, doi: 10.1029/98RS02381.
"""
function ioncollisionfrequency(z)
    return 4.54e9*exp(-0.15e-3*z)  # e-3 converts `z` to km
end

"""
    ioncollisionfrequency(z, cutoff_low)

When `z` is below `cutoff_low` (m), return 0.
"""
function ioncollisionfrequency(z, cutoff_low)
    if z > cutoff_low
        return ioncollisionfrequency(z)
    else
        return zero(promote_type(Float64, typeof(z)))
    end
end
