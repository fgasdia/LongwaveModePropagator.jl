export Constituent
export waitprofile, collisionprofile
export earthradius, speedoflight, fundamentalcharge, mₑ, μ₀, ϵ₀

const earthradius = 6369.427  # km, because it is used with heights
const speedoflight = 299792458  # m/s
const μ₀ = 4e-7π  # H/m
const ϵ₀ = 1/(μ₀*speedoflight^2)  # F/m
const mₑ = 9.10938356e-31  # kg
const fundamentalcharge = 1.602176621e-19  # C

struct Constituent{T<:Number, F<:Function, G<:Function}
    charge::T  # C
    mass::T  # kg
    numberdensity::F  # function that obtains number density (m⁻³)
    collisionfrequency::G  # function that obtains collision frequency
end

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

# struct Source{T} <: AbstractSource
#     freq::T
#     ω::T
#     k::T
#
#     function Source{T}(freq::T) where T<:Number
#         # f, ω, k
#         ω = 2π*freq
#         new(freq, ω, 1000ω/speedoflight) # NOTE: `k` is scalled to `km`
#     end
# end
# Source(freq::T) where T<:Number = Source{T}(freq)


struct Ground
    ϵᵣ::Float64
    σ::Float64
end

"""
Return the magnetic field parameters (`B`, `dcm`, `dcn`, `dcl`) at a specific location.
"""
function magneticfield()
    B = 0.0
    dcm = 0.0
    dcn = 0.0
    dcl = 0.0
    return B, dcm, dcn, dcl
end

"""
Calculate electron number density using Wait's profile (Wait 1964, Ferguson 1980). Returns
in m⁻³.

Note:
- `h′` must be in km
- `β` must be in km⁻¹
"""
function waitprofile(h, h′, β)
    return 1.43e13*exp(-0.15h′)*exp((β-0.15)*(h-h′))
end

"""
An exponential collision profile from Wait Spies 1964 (and MS 1976). Returns in collisions
per second.

Note:
- `h` must be in km
"""
function collisionprofile(h)
    return 1.86e11*exp(-0.15h)
end

"""
Based on prfl_exp.for and prfl_ennu.for.

NOT optimized for speed.
"""
function lwpccollisionprofile(h, h′, β)
    # Initialize (from prfl_exp.for)
    nrhtnu = 2
    htnu = Array{Float64}(undef, 51)
    algnu = Array{Float64}(undef, 51)

    htnu[2] = 0
    htnu[1] = 200
    algnu[2] = 25.925
    algnu[1] = -4.075

    ihtnu = 1
    if (h > htnu[ihtnu]) | (h <= htnu[ihtnu+1])
        if h > htnu[ihtnu]
            while (ihtnu > 1) & (h >= htnu[ihtnu-1])
                ihtnu -= 1
            end
            if ihtnu == 1
                # Extrapolate above top of profile
                j1 = 1
                j2 = 2
            else
                j1 = ihtnu - 1
                j2 = ihtnu
            end
        else
            while (ihtnu < nrhtnu) & (h <= htnu[ihtnu+1])
                ihtnu += 1
            end
            if ihtnu == nrhtnu
                # Extrapolate below the profile
                j1 = nrhtnu - 1
                j2 = nrhtnu
            else
                j1 = ihtnu
                j2 = ihtnu + 1
            end
        end
    else
        j1 = ihtnu
        j2 = ihtnu + 1
    end
    slope = (h - htnu[j1])/(htnu[j2]-htnu[j1])
    nu = exp(algnu[j1] + (algnu[j2]-algnu[j1])*slope)
end

"""
Based on prfl_exp.for and prfl_ennu.for.

NOT optimized at all.
"""
function lwpcdensityprofile(h, h′, β)
    # Initialize (from prfl_exp.for)
    nrhten = 2
    hten = Array{Float64}(undef, 51)
    algen = Array{Float64}(undef, 51)

    hten[2] = 0E
    hten[1] = 200
    algen[2] = 25.925-β*h′-9.4517306
    algen[1] = -4.075-β*h′-9.4517306+β*200

    ihten = 1
    if (h > hten[ihten]) | (h <= hten[ihten+1])
        if h > hten[ihten]
            while (ihten > 1) & (h >= hten[ihten-1])
                ihten -= 1
            end
            if ihten == 1
                # Extrapolate above top of profile
                j1 = 1
                j2 = 2
            else
                j1 = ihten - 1
                h2 = ihten
            end
        else
            while (ihten < nrhten) & (h <= hten[ihten+1])
                ihten += 1
            end
            if ihten == nrhten
                # Extrapolate below the profile
                j1 = nrhten - 1
                j2 = nrhten
            else
                j1 = ihten
                h2 = ihten + 1
            end
        end
    else
        j1 = ihten
        j2 = ihten + 1
    end
    slope = (h - hten[j1])/(hten[j2] - hten[j1])
    en = exp(algen[j1] + (algen[j2]-algen[j1])*slope)
end
