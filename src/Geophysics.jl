# CODATA 2018 NIST SP 961
const C0 = 299792458.0  # c, speed of light in vacuum, m/s
const U0 = 1.25663706212e-6  # μ₀, vacuum permeability, H/m
const E0 = 8.8541878128e-12  # ϵ₀, vacuum permittivity, F/m
const Z0 = 376.730313668  # Z₀, vacuum impedance, Ω

const QE = -1.602176634e-19  # qₑ, charge of an electron, C
const ME = 9.1093837015e-31  # mₑ, mass of an electron, kg

########

"""
    Species{F, G}

Ionosphere constituent `Species` with a `charge` (C), `mass` (kg),
`numberdensity(z)` (m⁻³), and `collisionfrequency(z)` (s⁻¹) where the latter two
should be callable functions of height `z`.
"""
struct Species{F, G}
    charge::Float64  # C
    mass::Float64  # kg
    numberdensity::F  # m⁻³
    collisionfrequency::G  # s⁻¹
end

########

"""
    BField

Background magnetic field vector of strength `B` (T) with direction cosines `dcl`, `dcm`,
and `dcn` corresponding to x, y, and z directions.
"""
struct BField
    B::Float64
    dcl::Float64
    dcm::Float64
    dcn::Float64

    function BField(B, dcl, dcm, dcn)
        if iszero(B)
            @warn "B field magnitude of exactly 0 is not supported. Setting B = 1e-15."
            B = 1e-15
        end
        new(B, dcl, dcm, dcn)
    end
end

"""
    BField(B, dip, azimuth)

Return a `BField` vector of strength `B` (T), `dip` angle (rad) from the
horizontal, and `azimuth` angle (rad) from the propagation direction, positive towards y.

Dip (inclination) is the angle made with the horizontal by the background magnetic field
lines. Positive dip corresponds to when B is directed downward into Earth. Therefore, it
ranges between +π in Earth's geographic northern hemisphere and -π in the southern
hemisphere.
"""
function BField(B, dip, azimuth)
    # BUG: How can we avoid this?
    # Look at the Booker quartic roots for dip angles -1,+1. They are way
    # outside the normal. The roots aren't accurately found, although I'm not
    # sure where the issue is.
    abs(dip) <= deg2rad(1) && @warn "magnetic dip angles between ±1° have known numerical issues."
    abs(dip) > π && @warn "magnetic dip angle should be in radians"
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

function isisotropic(b::BField)
    b.B <= 1e-12 && return true

    tolerance = deg2rad(0.15)

    bdip = dip(b)
    # BUG?: this is the opposite of what I'd first expect... it means the dip angle should be horizontal to be isotropic
    abs(bdip) < tolerance && return true

    baz = azimuth(b)
    (abs(baz - π/2) < tolerance || abs(baz - 3π/2) < tolerance) && return true

    return false
end

########
# Useful plasma functions

"""
    plasmafrequency(n, q, m)

Return the plasma frequency ``ωₚ = √(nq²/(ϵ₀m))`` for a ``cold'' plasma.
"""
plasmafrequency(n, q, m) = sqrt(n*q^2/(E0*m))
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

"Default ground conductivity indices from LWPC."
const GROUND = Dict(
    1=>Ground(5, 1e-5),
    2=>Ground(5, 3e-5),
    3=>Ground(10, 1e-4),
    4=>Ground(10, 3e-4),
    5=>Ground(15, 1e-3),
    6=>Ground(15, 3e-3),
    7=>Ground(15, 1e-2),
    8=>Ground(15, 3e-2),
    9=>Ground(15, 1e-1),
    10=>Ground(81, 4.0))

########

"""
    Fields

This `baremodule` allows scoped enum-like access to electric field components
`Ex`, `Ey`, and `Ez` as e.g., `Fields.Ex`.
"""
baremodule Fields
using Base: @enum
@enum Field Ex Ey Ez
end

########

"""
    waitprofile(z, h′, β; cutoff_low=0, threshold=1e12)

Return the electron number density in m⁻³ at altitude `z` (m) using Wait's
exponential profile with parameters `h′` (km) and `β` (km⁻¹).

``Nₑ = 1.43 × 10¹³ \\exp(-0.15 h') ⋅ \\exp[(β - 0.15)(z/1000 - h')]``

Optional Arguments:

    - When `z` is below `cutoff_low`, return 0.
    - When density is greater than `threshold`, return `threshold`

See also: [`electroncollisionfrequency`](@ref), [`ioncollisionfrequency`](@ref)

# References
[^Wait1964]: J. R. Wait and K. P. Spies, “Characteristics of the
    earth-ionosphere waveguide for VLF radio waves,” U.S. National Bureau of
    Standards, Boulder, CO, Technical Note 300, Dec. 1964.

[^Thomson1993]: N. R. Thomson, “Experimental daytime VLF ionospheric
    parameters,” Journal of Atmospheric and Terrestrial Physics, vol. 55, no. 2,
    pp. 173–184, Feb. 1993, doi: 10.1016/0021-9169(93)90122-F.

[^Cummer1998]: S. A. Cummer, U. S. Inan, and T. F. Bell, “Ionospheric D region
    remote sensing using VLF radio atmospherics,” Radio Science, vol. 33, no. 6,
    pp. 1781–1792, Nov. 1998, doi: 10.1029/98RS02381.
"""
function waitprofile(z, hp, β; cutoff_low=0, threshold=1e12)
    if z > cutoff_low
        # Using form with single `exp` for speed
        Ne = 1.43e13*exp(-0.15*hp - (β - 0.15)*(hp - z/1000))
        if Ne > threshold
            return oftype(Ne, threshold)
        else
            return Ne
        end
    else
        zero(promote_type(Float64, typeof(z)))
    end
end

"""
    electroncollisionfrequency(z)

Return the electron-neutral collision frequency (s⁻¹) at altitude `z` (m) based on Wait's
conductivity profile.

``νₑ(z) = 1.816 × 10¹¹ \\exp(-0.15e-3 z)``

See also: [`waitprofile`](@ref), [`ioncollisionfrequency`](@ref)

# References
[^Wait1964]: J. R. Wait and K. P. Spies, “Characteristics of the
earth-ionosphere waveguide for VLF radio waves,” U.S. National Bureau of
Standards, Boulder, CO, Technical Note 300, Dec. 1964.

[^Thomson1993]: N. R. Thomson, “Experimental daytime VLF ionospheric
parameters,” Journal of Atmospheric and Terrestrial Physics, vol. 55, no. 2,
pp. 173–184, Feb. 1993, doi: 10.1016/0021-9169(93)90122-F.

[^Cummer1998]: S. A. Cummer, U. S. Inan, and T. F. Bell, “Ionospheric D region
remote sensing using VLF radio atmospherics,” Radio Science, vol. 33, no. 6,
pp. 1781–1792, Nov. 1998, doi: 10.1029/98RS02381.
"""
function electroncollisionfrequency(z)
    return 1.816e11*exp(-0.15e-3*z)  # e-3 converts `z` to km
end

"""
    ioncollisionfrequency(z)

Return ion-neutral collision frequency (s⁻¹) at height `z` (m) from [^Morfitt1976].

``νᵢ(z) = 4.54 × 10⁹ \\exp(-0.15e-3 z)``

See also: [`waitprofile`](@ref), [`electroncollisionfrequency`](@ref)

# References
[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved
computer program for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere
waveguide,” Naval Electronics Laboratory Center, San Diego, CA, NELC/IR-77T,
Oct. 1976.

[^Cummer1998]: S. A. Cummer, U. S. Inan, and T. F. Bell, “Ionospheric D region
remote sensing using VLF radio atmospherics,” Radio Science, vol. 33, no. 6,
pp. 1781–1792, Nov. 1998, doi: 10.1029/98RS02381.
"""
function ioncollisionfrequency(z)
    return 4.54e9*exp(-0.15e-3*z)  # e-3 converts `z` to km
end
