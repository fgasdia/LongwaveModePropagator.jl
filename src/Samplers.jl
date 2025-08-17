"""
    enum Fields

Vector components of the electric field.

    - `Ex` is in the direction from the transmitter toward the receiver
    - `Ez` extends vertically upward into the ionosphere
    - `Ey` is perpendicular to the ``x-z`` plane in which the wavefields propagate,
        completing the right-handed coordinate system
    - `E` calculates `Ez`, `Ey`, and `Ex`
"""
@enumx Fields begin
    Ez
    Ey
    Ex
    E
end

function index(f::Fields.T)
    if f == Fields.Ez
        return 1
    elseif f == Fields.Ey
        return 2
    elseif f == Fields.Ex
        return 3
    elseif f == Fields.E
        return 1:3
    end
end
numcomponents(f::Fields.T) = length(index(f))

########

"""
    AbstractSampler

Abstract supertype for sampling fields in the waveguide.

Subtypes of AbstractSampler have a position in the guide and a field component.
"""
abstract type AbstractSampler{T} end

function distancechecks(d, stype; earthradius=LMPParams().earthradius)
    !issorted(d) && error("$stype distance must be sorted")
    minimum(d) < 0 && error("$stype distance must be positive")
    maximum(d) > (π*earthradius - 500e3) && @warn "antipode focusing effects not modeled"

    return true
end

"""
    Sampler{S} <: AbstractSampler{S}

`Sampler` types sample (measure) the electromagnetic field in the waveguide.

# Fields

- `distance::S`: ground distance from the transmitter in meters.
- `fieldcomponent::Fields.T`: field component measured by the `Sampler`.
- `altitude::Float64`: height above the ground in meters.
"""
struct Sampler{S} <: AbstractSampler{S}
    distance::S
    fieldcomponent::Fields.T
    altitude::Float64

    Sampler(d::S, fc, a) where S = distancechecks(d, Sampler) && new{S}(d, fc, a)
end
distance(s::Sampler) = s.distance
distance(s::Sampler,t::Transmitter) = s.distance
fieldcomponent(s::Sampler) = s.fieldcomponent
altitude(s::Sampler) = s.altitude

"""
    GroundSampler{S} <: AbstractSampler{S}

`GroundSamplers` are `Sampler` types with an altitude of zero.

# Fields

- `distance::S`: ground distance from the transmitter in meters.
- `fieldcomponent::Fields.T`: field component measured by the `GroundSampler`.
"""
struct GroundSampler{S} <: AbstractSampler{S}
    distance::S  # S is usually <: AbstractVector but could be a scalar
    fieldcomponent::Fields.T

    GroundSampler(d::S, fc) where S = distancechecks(d, GroundSampler) && new{S}(d, fc)
end
distance(g::GroundSampler) = g.distance
distance(g::GroundSampler,t::Transmitter) = g.distance
fieldcomponent(g::GroundSampler) = g.fieldcomponent
altitude(g::GroundSampler) = 0.0

"""
    Receiver{<:Antenna}

Represent a physical receiver.

A default `Receiver{VerticalDipole}` is returned with `Receiver()`.

# Fields

- `name::String = ""`: receiver name.
- `latitude::Float64 = 0.0`: geographic latitude in degrees.
- `longitude::Float64 = 0.0`: geographic longitude in degrees.
- `altitude::Float64 = 0.0`: receiver altitude in meters above the ground.
- `antenna::Antenna = VerticalDipole()`: receiver antenna.
"""
struct Receiver{A<:Antenna}
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
end

Receiver() = Receiver{VerticalDipole}("", 0.0, 0.0, 0.0, VerticalDipole())

altitude(r::Receiver) = r.altitude
distance(r::Receiver,t::Transmitter) = nothing   # TODO proj4 stuff
fieldcomponent(r::Receiver) = fieldcomponent(r.antenna)
