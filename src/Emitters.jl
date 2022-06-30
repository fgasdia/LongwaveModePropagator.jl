# Electromagnetic transmitters and related functions

wavenumber(f) = angular(f)/C0
wavelength(f) = C0/f
angular(f) = 2π*f

########

"""
    Emitter

Abstract type for types that energize the waveguide with electromagnetic energy.
"""
abstract type Emitter end

"""
    Transmitter{A<:Antenna} <: Emitter

Typical ground-based Transmitter.

# Fields

- `name::String`: transmitter name.
- `latitude::Float64`: transmitter geographic latitude in degrees.
- `longitude::Float64`: transmitter geographic longitude in degrees.
- `antenna::Antenna`: transmitter antenna.
- `frequency`: transmitter frequency in Hertz.
- `power::Float64`: transmit power in Watts.
"""
struct Transmitter{A<:Antenna} <: Emitter
    name::String
    latitude::Float64
    longitude::Float64
    antenna::A
    frequency::Float64
    power::Float64

    Transmitter{A}(n, la, lo, an, fr, po) where {A<:Antenna} = new(n, la, lo, an, fr, po)
end

altitude(t::Transmitter) = 0.0
frequency(t::Transmitter) = t.frequency
azimuth(t::Transmitter) = azimuth(t.antenna)
inclination(t::Transmitter) = inclination(t.antenna)
power(t::Transmitter) = t.power

"""
    Transmitter(name, lat, lon, frequency)

Return a `Transmitter` with a `VerticalDipole` antenna and transmit power of 1000 W.
"""
Transmitter(name, lat, lon, frequency) =
    Transmitter{VerticalDipole}(name, lat, lon, VerticalDipole(), frequency, 1000)

"""
    Transmitter(frequency)

Return a `Transmitter` with zeroed geographic position, 1000 W transmit power, and
`VerticalDipole` antenna.
"""
Transmitter(frequency) = Transmitter("", 0, 0, frequency)

"""
    Transmitter(antenna, frequency, power)

Return a `Transmitter` with zeroed geographic position.
"""
Transmitter(antenna::A, frequency, power) where {A<:Antenna} =
    Transmitter{A}("", 0, 0, antenna, frequency, power)
