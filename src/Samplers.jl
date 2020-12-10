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

`Sampler` types sample the field specified by `Fields.Field` at `altitude` and ground
`distance` from the transmitter.
"""
struct Sampler{T} <: AbstractSampler{T}
    distance::T
    altitude::Float64
    fieldcomponent::Fields.Field

    Sampler(d::T, fc) where T = distancechecks(d, Sampler) && new{T}(d, fc)
end
altitude(s::Sampler) = s.altitude
distance(s::Sampler) = s.distance
distance(s::Sampler,t::Transmitter) = s.distance
fieldcomponent(s::Sampler) = s.fieldcomponent

"""
    GroundSampler{T} <: AbstractSampler{T}

`GroundSamplers` are `Sampler` types with an altitude of zero.
"""
struct GroundSampler{T} <: AbstractSampler{T}
    distance::T  # T is usually <: AbstractVector but could be a scalar
    fieldcomponent::Fields.Field

    GroundSampler(d::T, fc) where T = distancechecks(d, GroundSampler) && new{T}(d, fc)
end
altitude(g::GroundSampler) = 0.0
distance(g::GroundSampler) = g.distance
distance(g::GroundSampler,t::Transmitter) = g.distance
fieldcomponent(g::GroundSampler) = g.fieldcomponent

"""
    Receiver{<:Antenna}

Represent a physical receiver.

# Fields

- `name::String`: receiver name.
- `latitude::Float64`: geographic latitude in degrees.
- `longitude::Float64`: geographic longitude in degrees.
- `altitude::Float64`: receiver altitude in meters above the ground.
- `antenna::Antenna`: receiver antenna.
"""
struct Receiver{A<:Antenna}
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
end

"""
    Receiver()

Return a default `Receiver{VerticalDipole}` with an empty name and zeroed geographic
position.
"""
Receiver() = Receiver{VerticalDipole}("", 0.0, 0.0, 0.0, VerticalDipole())

altitude(r::Receiver) = r.altitude
distance(r::Receiver,t::Transmitter) = nothing   # TODO proj4 stuff
fieldcomponent(r::Receiver) = fieldcomponent(r.antenna)
