export Constituent
export waitprofile, electroncollisionfrequency, ioncollisionfrequency
export earthradius, speedoflight, fundamentalcharge, electronmass, μ₀, ϵ₀

const earthradius = 6369427  # m
const speedoflight = 299792458  # m/s
const μ₀ = 4e-7π  # H/m
const ϵ₀ = 1/(μ₀*speedoflight^2)  # F/m
const electronmass = 9.10938356e-31  # kg
const fundamentalcharge = 1.602176621e-19  # C

struct Constituent{T<:Number, F<:Function, G<:Function}
    charge::T  # C
    mass::T  # kg
    numberdensity::F  # function that obtains number density (m⁻³)
    collisionfrequency::G  # function that obtains collision frequency
end

"""
`B` should be in Teslas
"""
struct BField{T}
    B::T
    dcl::T
    dcm::T
    dcn::T
end
function BField(B::T, dip::T, azimuth::T) where T<:Number
    k1 = cosd(dip)
    BField{T}(B, k1*cosd(azimuth), k1*sind(azimuth), -sind(dip))
end

struct Ground{S,T}
    ϵᵣ::S  # usually happens to be an int
    σ::T
end


"""
Calculate electron number density using Wait's profile (Wait 1964, Ferguson 1980). Returns
in m⁻³.

Note:
- `z` must be in m
- `h′` must be in km
- `β` must be in km⁻¹
"""
function waitprofile(z, h′, β)
    return 1.43e13*exp(-0.15h′)*exp((β - 0.15)*(z*1e-3 - h′))
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

"""
Exponential ion-neutral collision frequency from Morfitt and Shellman, 1976.
Returns in collisions per second. Also see Cummer et al 1998.

Note:
- `z` must be in m
"""
function ioncollisionfrequency(z)
    return 4.54e9*exp(-0.15e-3*z)  # converts `z` to km
end
