# TODO: Custom "range" of Frequencies to form a spectrum for, e.g. lightning sferic

"""
    Frequency

Electromagnetic wave defined by frequency `f`, but also carrying angular wave frequency `ω`,
wavenumber `k`, and wavelength `λ`, all in SI units.
"""
struct Frequency
    f::Float64
    ω::Float64
    k::Float64
    λ::Float64
end

"""
    Frequency(f)

Return `Frequency` given electromagnetic wave frequency `f` in Hertz.
"""
function Frequency(f)
    ω = 2π*f
    k = ω/C0
    λ = C0/f

    Frequency(f, ω, k, λ)
end
Frequency(f::Frequency) = f

########

"""
    Emitter

Abstract type for types that energize the waveguide with electromagnetic energy.
"""
abstract type Emitter end

function transmitterchecks(alt)
    alt < 0 && error("Transmitter altitude cannot be negative.")

    return true
end

"""
    Transmitter{A<:Antenna} <: Emitter

Typical ground-based Transmitter.

# Fields

- `name::String`: transmitter name.
- `latitude::Float64`: transmitter geographic latitude in degrees.
- `longitude::Float64`: transmitter geographic longitude in degrees.
- `antenna::Antenna`: transmitter antenna.
- `frequency::Frequency`: transmit frequency.
- `power::Float64`: transmit power in Watts.
"""
struct Transmitter{A<:Antenna} <: Emitter
    name::String
    latitude::Float64
    longitude::Float64
    antenna::A
    frequency::Frequency
    power::Float64

    Transmitter{A}(n, la, lo, an, fr, po) where {A<:Antenna} = new(n, la, lo, an, fr, po)
end

altitude(t::Transmitter) = 0.0
frequency(t::Transmitter) = t.frequency
azimuth(t::Transmitter) = azimuth(t.antenna)
inclination(t::Transmitter) = inclination(t.antenna)
power(t::Transmitter) = t.power

"""
    Transmitter(name, lat, lon, freq)

Return a `Transmitter` with a `VerticalDipole` antenna and transmit power of 1000 W.
"""
Transmitter(name, lat, lon, freq) =
    Transmitter{VerticalDipole}(name, lat, lon, VerticalDipole(), Frequency(freq), 1000)

"""
    Transmitter(freq)

Return a `Transmitter` with zeroed geographic position, 1000 W transmit power, and
`VerticalDipole` antenna.
"""
Transmitter(freq) = Transmitter("", 0, 0, freq)

"""
    Transmitter(freq, power, antenna)

Return a `Transmitter` with zeroed geographic position.
"""
Transmitter(antenna::A, freq::Frequency, power) where {A<:Antenna} =
    Transmitter{A}("", 0, 0, antenna, Frequency(freq), power)

"""
    AirborneTransmitter{A<:Antenna} <: Emitter

Similar to `Transmitter` except it also contains an `altitude` field in meters.
"""
struct AirborneTransmitter{A<:Antenna} <: Emitter
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
    frequency::Frequency
    power::Float64

    AirborneTransmitter{A}(n, la, lo, al, an, fr, po) where {A<:Antenna} =
        transmitterchecks(al) && new(n, la, lo, al, an, fr, po)
end

altitude(t::AirborneTransmitter) = t.altitude
frequency(t::AirborneTransmitter) = t.frequency
azimuth(t::AirborneTransmitter) = azimuth(t.antenna)
inclination(t::AirborneTransmitter) = inclination(t.antenna)
power(t::AirborneTransmitter) = t.power
