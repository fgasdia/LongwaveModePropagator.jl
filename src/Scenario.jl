export Receiver, Transmitter, Inputs, Source, ExcitationDipole

abstract type AbstractSource end

struct Receiver
    name::String
    lat::Float64
    lon::Float64
    alt::Float64
end

@enum ExcitationDipole begin
    vertical
    horizontal_endon
    horizontal_broadside
end

struct Source{T} <: AbstractSource
    freq::T
    ω::T
    k::T
    λ::T
    exciter::ExcitationDipole

    function Source{T}(freq::T, exciter::ExcitationDipole) where T<:Real
        # f, ω, k
        ω = 2π*freq
        k = ω/speedoflight
        λ = speedoflight/freq

        new(freq, ω, k, λ, exciter)
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
