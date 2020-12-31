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

"""
    AbstractSampler

Abstract supertype for sampling fields in the waveguide.

Subtypes of AbstractSampler have a position in the guide and a `Fields.Field`.
"""
abstract type AbstractSampler{T} end

function distancechecks(d, stype; earthradius=LMPParams().earthradius)
    !issorted(d) && error("$stype distance must be sorted")
    minimum(d) < 0 && error("$stype distance must be positive")
    maximum(d) > (π*earthradius - 500e3) && @warn "antipode focusing effects not modeled"

    return true
end

"""
    Sampler{T} <: AbstractSampler{T}

`Sampler` types sample (measure) the electromagnetic field in the waveguide.

# Fields

- `distance::T`: ground distance from the transmitter in meters.
- `fieldcomponent::Fields.Field`: field component measured by the `Sampler`.
- `altitude::Float64`: height above the ground in meters.
"""
struct Sampler{T} <: AbstractSampler{T}
    distance::T
    fieldcomponent::Fields.Field
    altitude::Float64

    Sampler(d::T, fc, a) where T = distancechecks(d, Sampler) && new{T}(d, fc, a)
end
distance(s::Sampler) = s.distance
distance(s::Sampler,t::Transmitter) = s.distance
fieldcomponent(s::Sampler) = s.fieldcomponent
altitude(s::Sampler) = s.altitude

"""
    GroundSampler{T} <: AbstractSampler{T}

`GroundSamplers` are `Sampler` types with an altitude of zero.

# Fields

- `distance::T`: ground distance from the transmitter in meters.
- `fieldcomponent::Fields.Field`: field component measured by the `GroundSampler`.
"""
struct GroundSampler{T} <: AbstractSampler{T}
    distance::T  # T is usually <: AbstractVector but could be a scalar
    fieldcomponent::Fields.Field

    GroundSampler(d::T, fc) where T = distancechecks(d, GroundSampler) && new{T}(d, fc)
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
