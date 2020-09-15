"""
    AbstractSampler

Abstract supertype for sampling fields in the waveguide.

Subtypes of AbstractSampler should have a position in the guide and a FieldComponent
"""
abstract type AbstractSampler end

function distancechecks(d, stype)
    !issorted(d) && error("$t distance must be sorted")
    minimum(d) < 0 && error("$t distance must be positive")
    maximum(d) > (π*EARTHRADIUS - 500e3) && @warn "antipode focusing effects not modeled"

    return true
end

struct Sampler{R<:AbstractVector} <: AbstractSampler
    distance::R
    altitude::Float64
    fieldcomponent::FieldComponent

    Sampler{R}(d, fc) where {R<:AbstractVector} = distancechecks(d, Sampler) && new(d, fc)
end
altitude(s::Sampler) = s.altitude
distance(s::Sampler) = s.distance
distance(s::Sampler,t::Transmitter) = s.distance
fieldcomponent(s::Sampler) = s.fieldcomponent

struct GroundSampler{R<:AbstractVector} <: AbstractSampler
    distance::R
    fieldcomponent::FieldComponent

    GroundSampler{R}(d, fc) where {R<:AbstractVector} = distancechecks(d, GroundSampler) && new(d, fc)
end
GroundSampler(d, fc) = GroundSampler{typeof(d)}(d, fc)
altitude(g::GroundSampler) = 0
distance(g::GroundSampler) = g.distance
distance(g::GroundSampler,t::Transmitter) = g.distance
fieldcomponent(g::GroundSampler) = g.fieldcomponent

struct Receiver{A<:Antenna} <: AbstractSampler
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
end
altitude(r::Receiver) = r.altitude
distance(r::Receiver,t::Transmitter) = nothing   # TODO proj4 stuff
fieldcomponent(r::Receiver) = fieldcomponent(r.antenna)
