# TODO: Custom "range" of Frequencies to form a spectrum for, e.g. lightning sferic

"""
    Frequency

Electromagnetic wave defined by frequency `f`, but also carrying angular wave
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
    AirborneTransmitter{A<:Antenna} <: Emitter

Similar to `Transmitter` except it also contains an `altitude` field.
"""
struct AirborneTransmitter{A<:Antenna} <: Emitter
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
    frequency::Frequency
    power::Float64

    AirborneTransmitter{A}(n, la, lo, al, an, fr, po) where {A<:Antenna} = transmitterchecks(al) && new(n, la, lo, al, an, fr, po)
end


"""
    Transmitter{A<:Antenna} <: Emitter

Typical ground-based Transmitter with a `name`, `latitude`, `longitude`,
`antenna`, `frequency`, and transmit `power`.

`altitude(Transmitter)` returns ``0.0``.
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

"""
    Transmitter(freq)

Generate a default transmitter with zeroed geographic position, 1000 W transmit
power, `VerticalDipole` antenna, at frequency `freq`.
"""
Transmitter(freq) = Transmitter("", 0, 0, freq)

"""
    Transmitter(name, lat, lon, freq)

Generate a default transmitter with a `name`, `lat`, `lon`, and frequency `freq`.

This will have a `VerticalDipole` antenna and transmit power of 1000 W.
"""
Transmitter(name, lat, lon, freq) =
    Transmitter{VerticalDipole}(name, lat, lon, VerticalDipole(), Frequency(freq), 1000)

"""
    Transmitter(freq, power, antenna)

Generate a default transmitter with frequency `freq`, transmit `power`, and
`antenna` and zeroed geographic position.
"""
Transmitter(freq, power, antenna::A) where {A<:Antenna} =
    Transmitter{A}("", 0, 0, antenna, Frequency(freq), power)

altitude(t::Transmitter) = 0.0
frequency(t::Transmitter) = t.frequency
azimuth(t::Transmitter) = azimuth(t.antenna)
inclination(t::Transmitter) = inclination(t.antenna)
power(t::Transmitter) = t.power
