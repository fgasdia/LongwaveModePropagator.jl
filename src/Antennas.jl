"""
    Antenna

Abstract supertype for all antenna-like types. 
"""
abstract type Antenna end

"""
    AbstractDipole <: Antenna

Abstract type for dipole antennas.

Subtypes of `AbstractDipole` have an orientation `azimuth_angle` (``ϕ``) and
`inclination_angle` (``γ``) where

| [rad] |     ϕ     |     γ      |
|:-----:|:---------:|:----------:|
|   0   | end fire  | vertical   |
|  π/2  | broadside | horizontal |

!!! note

    Angles describe the orientation of the antenna, not the radiation pattern.
"""
abstract type AbstractDipole <: Antenna end

"""
    azimuth(d::AbstractDipole)

Return the azimuth angle ``ϕ`` of `d` measured towards positive `y` direction from `x`.

See also: [`AbstractDipole`](@ref)
"""
azimuth

"""
    inclination(d::AbstractDipole)

Return the inclination (dip) angle ``γ`` of `d` measured from vertical.

See also: [`AbstractDipole`](@ref)
"""
inclination

"""
    Dipole

Dipole antenna with arbitrary orientation described by `azimuth_angle` and
`inclination_angle` from vertical.

See also: [`AbstractDipole`](@ref)
"""
struct Dipole <: AbstractDipole
    azimuth_angle::Float64
    inclination_angle::Float64
end
azimuth(d::Dipole) = d.azimuth_angle
inclination(d::Dipole) = d.inclination_angle

"""
    VerticalDipole

Dipole antenna with inclination angle ``γ = 0`` from the vertical.
"""
struct VerticalDipole <: AbstractDipole end
azimuth(d::VerticalDipole) = 0.0
inclination(d::VerticalDipole) = 0.0

"""
    HorizontalDipole

Dipole antenna with inclination angle ``γ = π/2`` from the vertical and `azimuth_angle`
orientation ``ϕ`` where:

|  ϕ [rad] | orientation |
|:--------:|:-----------:|
|     0    |  end fire   |
|    π/2   |  broadside  |

"""
struct HorizontalDipole <: AbstractDipole
    azimuth_angle::Float64  # radians
end
azimuth(d::HorizontalDipole) = d.azimuth_angle
inclination(d::HorizontalDipole) = π/2
