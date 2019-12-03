export Receiver, Transmitter, Inputs, Source, ExcitationDipole

abstract type Exciter end
altitude(exc::Exciter) = exc.alt

"""
    VerticalDipole

Dipole antenna emitter at `alt` altitude in meters above ground.
"""
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
        k = ω/c₀
        λ = c₀/freq

        new{T}(freq, ω, k, λ, exciter)
    end
end
Source(freq) = Source(float(freq), VerticalDipole())
Source(freq, exciter::Exciter) = Source(float(freq), exciter)
source(s::Source) = s

struct Transmitter{S,T} <: AbstractSource
    name::String
    lat::S
    lon::S
    source::Source{T}
end
Transmitter(name, lat, lon, alt) = Transmitter(name, lat, lon, VerticalDipole(alt))
source(xmtr::Transmitter) = xmtr.source

abstract type AbstractReceiver end
altitude(r::AbstractReceiver) = r.alt

struct Receiver{S<:Real,T<:Real} <: AbstractReceiver
    name::String
    lat::S
    lon::S
    alt::T
end
