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

Ionosphere constituent `Species`.

# Fields

- `charge::Float64`: signed species charged in Coulombs.
- `mass::Float64`: species mass in kilograms.
- `numberdensity::F`: a callable that returns number density in number per cubic meter as a
    function of height in meters.
- `collisionfrequency::G`: a callable that returns the collision frequency in collisions per
    second as a function of height in meters.
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

Background magnetic field vector of strength `B` in Tesla with direction cosines `dcl`,
`dcm`, and `dcn` corresponding to ``x``, ``y``, and ``z`` directions parallel,
perpendicular, and up to the waveguide.
"""
struct BField
    B::Float64
    dcl::Float64
    dcm::Float64
    dcn::Float64

    function BField(B, dcl, dcm, dcn)
        if iszero(B)
            @warn "BField magnitude of exactly 0 is not supported. Setting B = 1e-16."
            B = 1e-16
        end
        new(B, dcl, dcm, dcn)
    end
end

"""
    BField(B, dip, azimuth)

Return a `BField` of strength `B` in Tesla, `dip` angle in radians from the horizontal, and
`azimuth` angle in radians from the propagation direction, positive towards ``y``.

Dip (inclination) is the angle made with the horizontal by the background magnetic field
lines. Positive dip corresponds to when the magnetic field vector is directed downward into
Earth. Therefore, it ranges between ``+π`` in Earth's geographic northern hemisphere and
``-π`` in the southern hemisphere.
"""
function BField(B, dip, azimuth)
    # BUG: How can we avoid this?
    # Look at the Booker quartic roots for dip angles -1,+1. They are way
    # outside the normal. The roots aren't accurately found, although I'm not
    # sure where the issue is.
    abs(dip) <= deg2rad(1) &&
        @warn "magnetic dip angles between ±1° have known numerical issues."
    abs(dip) > π &&
        @warn "magnetic dip angle should be in radians."
    abs(azimuth) > 2π &&
        @warn "magnetic azimuth angle should be in radians."

    Sdip, Cdip = sincos(dip)
    Saz, Caz = sincos(azimuth)
    BField(B, Cdip*Caz, Cdip*Saz, -Sdip)
end

"""
    dip(b::BField)

Return the dip angle in radians from `b`.
"""
dip(b::BField) = -asin(b.dcn)

"""
    azimuth(b::BField)

Return the azimuth angle in radians from `b`.
"""
azimuth(b::BField) = atan(b.dcm,b.dcl)

function isisotropic(b::BField)
    b.B <= 1e-12 && return true

    tolerance = deg2rad(0.15)

    bdip = dip(b)
    abs(bdip) < tolerance && return true

    baz = azimuth(b)
    (abs(baz - π/2) < tolerance || abs(baz - 3π/2) < tolerance) && return true

    return false
end

########

"""
    Ground(ϵᵣ, σ)

Isotropic earth ground characterized by a relative permittivity `ϵᵣ` and a conductivity `σ`
in Siemens per meter.
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

This `baremodule` allows scoped enum-like access to electric field components `Ex`, `Ey`,
and `Ez`.

# Examples

```jldoctest; setup=:(using LongwaveModePropagator.Fields)
julia> Fields.Ex
Ex::Field = 0
julia> Fields.Ey
Ey::Field = 1
```
"""
baremodule Fields
using Base: @enum
@enum Field Ex Ey Ez
end

########

@doc raw"""
    waitprofile(z, h′, β; cutoff_low=0, threshold=1e12)

Compute the electron number density in electrons per cubic meter at altitude `z` in meters
using Wait's exponential profile [[Wait1964](@cite); [Thomson1993](@cite)] with parameters
`h′` in kilometers and `β` in inverse kilometers.

The profile is:
```math
Nₑ = 1.43 × 10¹³ \exp(-0.15 h') \exp[(β - 0.15)(z/1000 - h')]
```

Optional arguments:
- `cutoff_low=0`: when `z` is below `cutoff_low`, return zero.
- `threshold=1e12`: when density is greater than `threshold`, return `threshold`.

See also: [`electroncollisionfrequency`](@ref), [`ioncollisionfrequency`](@ref)

# References

[Wait1964]: J. R. Wait and K. P. Spies, “Characteristics of the earth-ionosphere waveguide
    for VLF radio waves,” U.S. National Bureau of Standards, Boulder, CO, Technical Note
    300, Dec. 1964.

[Thomson1993]: N. R. Thomson, “Experimental daytime VLF ionospheric parameters,” Journal of
    Atmospheric and Terrestrial Physics, vol. 55, no. 2, pp. 173–184, Feb. 1993,
    doi: 10.1016/0021-9169(93)90122-F.
"""
function waitprofile(z, hp, β; cutoff_low=0, threshold=1e12)
    if z > cutoff_low
        # Using form with single `exp` for speed
        Ne = 1.43e13*exp(-0.15*hp - (β - 0.15)*(hp - z*0.001))  # converting `z` to km
        if Ne > threshold
            return oftype(Ne, threshold)
        else
            return Ne
        end
    else
        zero(promote_type(Float64, typeof(z)))
    end
end

@doc raw"""
    electroncollisionfrequency(z)

Compute the electron-neutral collision frequency in collisions per second at height `z` in
meters based on Wait's conductivity profile [[Wait1964](@cite); [Thomson1993](@cite)].

The profile is:
```math
νₑ(z) = 1.816 × 10¹¹ \exp(-0.15⋅z/1000)
```

See also: [`waitprofile`](@ref), [`ioncollisionfrequency`](@ref)

# References

[Wait1964]: J. R. Wait and K. P. Spies, “Characteristics of the earth-ionosphere waveguide
    for VLF radio waves,” U.S. National Bureau of Standards, Boulder, CO, Technical Note
    300, Dec. 1964.

[Thomson1993]: N. R. Thomson, “Experimental daytime VLF ionospheric parameters,” Journal of
    Atmospheric and Terrestrial Physics, vol. 55, no. 2, pp. 173–184, Feb. 1993,
    doi: 10.1016/0021-9169(93)90122-F.
"""
function electroncollisionfrequency(z)
    return 1.816e11*exp(-0.15e-3*z)  # e-3 converts `z` to km
end

@doc raw"""
    ioncollisionfrequency(z)

Compute the ion-neutral collision frequency in collisions per second at height `z` in meters
from [[Morfitt1976](@cite)].

The profile is:
```math
νᵢ(z) = 4.54 × 10⁹ \exp(-0.15⋅z/1000)
```

See also: [`waitprofile`](@ref), [`electroncollisionfrequency`](@ref)

# References

[Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function ioncollisionfrequency(z)
    return 4.54e9*exp(-0.15e-3*z)  # e-3 converts `z` to km
end

########
# Useful plasma functions

"""
    plasmafrequency(n, q, m)
    plasmafrequency(z, s::Species)

Return the plasma frequency ``ωₚ = √(nq²/(ϵ₀m))`` for a ``cold'' plasma.
"""
plasmafrequency(n, q, m) = sqrt(n*q^2/(E0*m))
plasmafrequency(z, s::Species) = plasmafrequency(s.numberdensity(z), s.charge, s.mass)

"""
    gyrofrequency(q, m, B)
    gyrofrequency(s::Species, b::BField)

Return the _signed_ gyrofrequency ``Ω = qB/m``.
"""
gyrofrequency(q, m, B) = q*B/m
gyrofrequency(s::Species, b::BField) = gyrofrequency(s.charge, s.mass, b.B)
