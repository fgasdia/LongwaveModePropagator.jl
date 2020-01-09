export Receiver, Transmitter, Inputs, Source, ExcitationDipole

"""
    Exciter

Abstract type for types that energize the waveguide with electromagnetic energy.
"""
abstract type Exciter end
# exciter(exc::Exciter) = exc

"""
    Antenna

Abstract type for antenna designs.

Subtypes of `Antenna` have an orientation `azimuth_angle` (``ϕ``) and
`elevation_angle` (``γ``) where

|     |     ϕ     |     γ      |
| --- | --------- | ---------- |
|  0  | end fire  | vertical   |
| π/2 | broadside | horizontal |

"""
abstract type Antenna end

abstract type AbstractDipole <: Antenna end

"""
    Dipole

Dipole antenna with arbitrary `azimuth_angle` and `elevation_angle` orientation.
"""
struct Dipole <: AbstractDipole
    azimuth_angle::Float64
    elevation_angle::Float64
end
azimuth(d::Dipole) = d.azimuth_angle
elevation(d::Dipole) = d.elevation_angle
fieldcomponent(d::Dipole) = FC_ALL

"""
    VerticalDipole

Dipole antenna with elevation angle ``γ = 0``.
"""
struct VerticalDipole <: AbstractDipole end
azimuth(d::VerticalDipole) = 0
elevation(d::VerticalDipole) = 0
fieldcomponent(d::VerticalDipole) = FC_Ez

"""
    HorizontalDipole

Dipole antenna with elevation angle ``γ = π/2`` and `azimuth_angle` orientation ``ϕ`` where

|     |     ϕ     |
| --- | --------- |
|  0  | end fire  |
| π/2 | broadside |
"""
struct HorizontalDipole <: AbstractDipole
    azimuth_angle::Float64  # radians
end
azimuth(d::HorizontalDipole) = d.azimuth_angle
elevation(d::HorizontalDipole) = π/2
fieldcomponent(d::HorizontalDipole) = FC_ALL

"""
    Frequency

Electromagnetic wave component defined by frequency `f`, but also carrying angular wave
frequency, wavenumber, and wavelength.
"""
struct Frequency
    f::Float64
    ω::Float64
    k::Float64
    λ::Float64
end
function Frequency(f)
    ω = 2π*f
    k = ω/c₀
    λ = c₀/f

    Frequency(f, ω, k, λ)
end

# TODO: Custom "range" of Frequencies to form a spectrum for, e.g. lightning sferic

"""
    Transmitter
"""
struct Transmitter{A<:Antenna} <: Exciter
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
    frequency::Frequency
    power::Float64
end
Transmitter(name, lat, lon, freq) = Transmitter(name, lat, lon, 0, VerticalDipole(),
                                                Frequency(freq), 1)
altitude(t::Transmitter) = t.altitude
frequency(t::Transmitter) = t.frequency
azimuth(t::Transmitter) = azimuth(t.antenna)
elevation(t::Transmitter) = elevation(t.antenna)

"""
    AbstractSampler

Abstract supertype for sampling fields in the waveguide.

Subtypes of AbstractSampler should have a position in the guide.
"""
abstract type AbstractSampler end

struct Receiver{A<:Antenna} <: AbstractSampler
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
end
altitude(r::Receiver) = r.altitude
