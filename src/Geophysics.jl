export Constituent
export waitprofile, electroncollisionfrequency, ioncollisionfrequency
export Rₑ, c₀, μ₀, ϵ₀, H

# CODATA 2018 NIST SP 961
const Rₑ = 6366e3  # earth radius, m
const c₀ = 299792458  # speed of light in vacuum, m/s
const μ₀ = 1.25663706212e-6  # vacuum permeability, H/m
const ϵ₀ = 1/(μ₀*c₀^2)  # vacuum permittivity, F/m
const Z₀ = 376.730313668  # vacuum impedance, Ω
const H = 50e3  # reference height for Earth curvature, m
# const H = 0

struct Constituent{F<:Function, G<:Function}
    charge::Float64  # C
    mass::Float64  # kg
    numberdensity::F  # function that obtains number density (m⁻³)
    collisionfrequency::G  # function that obtains collision frequency
end

"""
`B` should be in Teslas, angles in radians!
"""
struct BField
    B::Float64
    dcl::Float64
    dcm::Float64
    dcn::Float64
end
function BField(B, dip, azimuth)
    Sdip, Cdip = sincos(dip)
    Saz, Caz = sincos(azimuth)
    BField(B, Cdip*Caz, Cdip*Saz, -Sdip)
end
# TODO: Function that returns dip and azimuth from direction cosines

struct Ground
    ϵᵣ::Int
    σ::Float64
end

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

struct EigenAngle{T}
    θ::T  # radians, because df/dθ are in radians
    cosθ::T
    sinθ::T
    cos²θ::T
    sin²θ::T

    function EigenAngle{T}(θ::T) where T <: Number
        rθ, iθ = reim(θ)
        (abs(rθ) > 2π) || (abs(iθ) > 2π) && @warn "θ in radians?"
        S, C = sincos(θ)
        C² = C^2
        S² = 1 - C²
        new(θ, C, S, C², S²)
    end
end
EigenAngle(θ::T) where T <: Number = EigenAngle{T}(θ)
Base.eltype(::Type{EigenAngle{T}}) where {T} = T


"""
Calculate electron number density using Wait's profile (Wait 1964, Ferguson 1980). Returns
in m⁻³.

``Nₑ = 1.43e13*exp(-0.15*h′)*exp((β - 0.15)*(z/1000 - h′))``

Note:
- `z` must be in m
- `h′` must be in km
- `β` must be in km⁻¹
"""
function waitprofile(z, h′, β)
    # Using form with single `exp` for speed
    return 1.43e13*exp(-0.15h′ - (β-0.15)*(h′ - z/1000))
end
function waitprofile(z, h′, β, cutoff_low)
    if z > cutoff_low
        return waitprofile(z, h′, β)
    else
        return zero(Float64)
    end
end

"""
Exponential electron-neutral collision frequency from Wait Spies 1964.
Returns in collisions per second. Also see Cummer et al 1998.

Note:
- `z` must be in m
"""
function electroncollisionfrequency(z)
    return 1.86e11*exp(-0.15e-3*z)  # converts `z` to km
end
function electroncollisionfrequency(z, cutoff_low)
    if z > cutoff_low
        return electroncollisionfrequency(z)
    else
        return zero(Float64)
    end
end

"""
Exponential ion-neutral collision frequency from Morfitt and Shellman, 1976.
Returns in collisions per second. Also see Cummer et al 1998.

Note:
- `z` must be in m
"""
function ioncollisionfrequency(z)
    return 4.54e9*exp(-0.15e-3*z)  # converts `z` to km
end
