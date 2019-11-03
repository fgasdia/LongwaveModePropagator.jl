export Receiver, Transmitter, Inputs, Source, ExcitationDipole

struct Receiver
    name::String
    lat::Float64
    lon::Float64
    alt::Float64
end

abstract type Exciter end

struct VerticalDipole <: Exciter
    alt::Float64  # m
end
VerticalDipole() = VerticalDipole(0.0)

"""
    HorizontalDipole

`ψ` is the angle between the direction of the horizontal dipole and the direction of
propagation, i.e. ``ψ=0`` represents end fire while ``ψ=π/2`` represents broadside launching.
"""
struct HorizontalDipole <: Exciter
    ψ::Float64  # radians
    alt::Float64  # m
end
HorizontalDipole(ψ) = HorizontalDipole(ψ, 0.0)

abstract type AbstractSource end

exciter(s::AbstractSource) = source(s).exciter
k(s::AbstractSource) = source(s).k

struct Source{T<:Number} <: AbstractSource
    freq::T
    ω::T
    k::T
    λ::T
    exciter::Exciter

    function Source(freq::T, exciter::Exciter) where T <: Number
        ω = 2π*freq
        k = ω/speedoflight
        λ = speedoflight/freq

        new{T}(freq, ω, k, λ, exciter)
    end
end
Source(freq) = Source(freq, VerticalDipole())
Source(freq, exciter::Exciter) = Source(freq, exciter)
source(s::Source) = s

struct Transmitter{S,T} <: AbstractSource
    name::String
    lat::S
    lon::S
    source::Source{T}
end
Transmitter(name, lat, lon, alt) = Transmitter(name, lat, lon, VerticalDipole(alt))
source(xmtr::Transmitter) = xmtr.source
