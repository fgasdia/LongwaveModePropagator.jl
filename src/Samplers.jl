"""
    AbstractSampler

Abstract supertype for sampling fields in the waveguide.

Subtypes of AbstractSampler have a position in the guide and a `FieldComponent`.
"""
abstract type AbstractSampler end

function distancechecks(d, stype)
    !issorted(d) && error("$stype distance must be sorted")
    minimum(d) < 0 && error("$stype distance must be positive")
    maximum(d) > (π*EARTHRADIUS - 500e3) && @warn "antipode focusing effects not modeled"

    return true
end

"""
    Sampler{R} <: AbstractSampler

`Sampler` types sample the field specified by `fieldcomponent` at `altitude` and
distance(s) `distance` from the transmitter.
"""
struct Sampler{R} <: AbstractSampler
    distance::R
    altitude::Float64
    fieldcomponent::FieldComponent

    Sampler(d::R, fc) where R = distancechecks(d, Sampler) && new{R}(d, fc)
end
altitude(s::Sampler) = s.altitude
distance(s::Sampler) = s.distance
distance(s::Sampler,t::Transmitter) = s.distance
fieldcomponent(s::Sampler) = s.fieldcomponent

"""
    GroundSampler{R} <: AbstractSampler

Ground samplers are designed to sample the field specified by `fieldcomponent`
at distance(s) `distance` along the ground from the transmitter.
"""
struct GroundSampler{R} <: AbstractSampler
    distance::R  # R is usually <: AbstractVector but could be a scalar
    fieldcomponent::FieldComponent

    GroundSampler(d::R, fc) where R = distancechecks(d, GroundSampler) && new{R}(d, fc)
end
altitude(g::GroundSampler) = 0.0
distance(g::GroundSampler) = g.distance
distance(g::GroundSampler,t::Transmitter) = g.distance
fieldcomponent(g::GroundSampler) = g.fieldcomponent

"""
    Receiver{<:Antenna} <: AbstractSampler

Represent a physical receiver with a `name`, `latitude`, `longitude`, `altitude`,
and `antenna`.
"""
struct Receiver{A<:Antenna} <: AbstractSampler
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
end

"""
    Receiver()

Returns a default `Receiver{VerticalDipole}` with an empty name and zeroed
geographic position.
"""
Receiver() = Receiver{VerticalDipole}("", 0.0, 0.0, 0.0, VerticalDipole())

altitude(r::Receiver) = r.altitude
distance(r::Receiver,t::Transmitter) = nothing   # TODO proj4 stuff
fieldcomponent(r::Receiver) = fieldcomponent(r.antenna)
