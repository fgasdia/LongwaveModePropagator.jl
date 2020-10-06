"""
    Antenna

Abstract type for antenna designs.
"""
abstract type Antenna end

"""
    AbstractDipole <: Antenna

Abstract type for dipole antenna designs.

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
    Dipole

Dipole antenna with arbitrary orientation described by `azimuth_angle` and
`inclination_angle` from vertical.
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

`azimuth(VerticalDipole())` returns ``0.0`` but physically this has
little meaning.
"""
struct VerticalDipole <: AbstractDipole end
azimuth(d::VerticalDipole) = 0.0
inclination(d::VerticalDipole) = 0.0

"""
    HorizontalDipole

Dipole antenna with inclination angle ``γ = π/2`` from the vertical and
`azimuth_angle` orientation ``ϕ`` where

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
