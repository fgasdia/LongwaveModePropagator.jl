"""
    AbstractSampler

Abstract supertype for sampling fields in the waveguide.

Subtypes of AbstractSampler should have a position in the guide and a FieldComponent
"""
abstract type AbstractSampler end

struct GroundSampler{R<:AbstractArray} <: AbstractSampler
    distance::R
    fieldcomponent::FieldComponent
end
altitude(g::GroundSampler) = 0
distance(g::GroundSampler) = g.distance
fieldcomponent(g::GroundSampler) = g.fieldcomponent

struct Receiver{A<:Antenna} <: AbstractSampler
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    antenna::A
end
altitude(r::Receiver) = r.altitude
fieldcomponent(r::Receiver) = fieldcomponent(r.antenna)

"""
Calculate distance from emitter
how do I do this without requiring to pass the transmitter?
may need a wrapper or a distance "flag" and if statement to see if distance needs to be
calculated

i.e. if isnothing(distance(r)) ... else
"""
distance(r::Receiver) = nothing
