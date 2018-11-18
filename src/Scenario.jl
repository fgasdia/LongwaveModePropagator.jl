export Receiver, Transmitter, Inputs, Source

abstract type AbstractSource end

struct Receiver
    name::String
    lat::Float64
    lon::Float64
    alt::Float64
end

struct Source{T} <: AbstractSource
    freq::T
    ω::T
    k::T

    function Source{T}(freq::T) where T<:Number
        # f, ω, k
        ω = 2π*freq
        new(freq, ω, ω/speedoflight)
    end
end
Source(freq::T) where T<:Number = Source{T}(freq)

struct Transmitter{T} <: AbstractSource
    name::String
    lat::T
    lon::T
    alt::T
    ω::T
    freq::T
    k::T
end

# TEMP: Only need to be mutable and have internal constructor for incomplete initialization
mutable struct Inputs
    transmitter::Transmitter
    receiver::Receiver
    maxrange::Float64
    deltarange::Float64
    topheight::Float64
    bottomheight::Float64
    Inputs() = new()
end
