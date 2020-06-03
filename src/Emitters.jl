# TODO: Custom "range" of Frequencies to form a spectrum for, e.g. lightning sferic

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
Frequency(f::Frequency) = f

########

"""
    Emitter

Abstract type for types that energize the waveguide with electromagnetic energy.
"""
abstract type Emitter end
# exciter(exc::Exciter) = exc

function transmitterchecks(alt)
    alt < 0 && error("Transmitter altitude cannot be negative.")

    return true
end

"""
    Transmitter
"""
struct Transmitter{A<:Antenna} <: Emitter
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
    frequency::Frequency
    power::Float64

    Transmitter{A}(n, la, lo, al, an, fr, po) where {A<:Antenna} = transmitterchecks(al) && new(n, la, lo, al, an, fr, po)
end
Transmitter(freq) = Transmitter("", 0, 0, freq)
Transmitter(name, lat, lon, freq) =
    Transmitter{VerticalDipole}(name, lat, lon, 0, VerticalDipole(), Frequency(freq), 1000)
Transmitter(freq, power, antenna::A) where {A<:Antenna} =
    Transmitter{A}("", 0, 0, 0, antenna, Frequency(freq), power)

altitude(t::Transmitter) = t.altitude
frequency(t::Transmitter) = t.frequency
azimuth(t::Transmitter) = azimuth(t.antenna)
elevation(t::Transmitter) = elevation(t.antenna)
power(t::Transmitter) = t.power
