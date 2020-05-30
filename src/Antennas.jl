"""
    Antenna

Abstract type for antenna designs.

Subtypes of `Antenna` have an orientation `azimuth_angle` (``ϕ``) and
`elevation_angle` (``γ``) where

|     |     ϕ     |     γ      |
| --- | --------- | ---------- |
|  0  | end fire  | vertical   |
| π/2 | broadside | horizontal |

"""
abstract type Antenna end

abstract type AbstractDipole <: Antenna end

"""
    Dipole

Dipole antenna with arbitrary `azimuth_angle` and `elevation_angle` orientation.

Note: these are physical orientation angles of the antenna, not the radiation pattern.
"""
struct Dipole <: AbstractDipole
    azimuth_angle::Float64
    elevation_angle::Float64
end
azimuth(d::Dipole) = d.azimuth_angle
elevation(d::Dipole) = d.elevation_angle

"""
    VerticalDipole

Dipole antenna with elevation angle ``γ = 0``.
"""
struct VerticalDipole <: AbstractDipole end
azimuth(d::VerticalDipole) = 0
elevation(d::VerticalDipole) = π/2

"""
    HorizontalDipole

Dipole antenna with elevation angle ``γ = π/2`` and `azimuth_angle` orientation ``ϕ`` where

|     |     ϕ     |
| --- | --------- |
|  0  | end fire  |
| π/2 | broadside |
"""
struct HorizontalDipole <: AbstractDipole
    azimuth_angle::Float64  # radians
end
azimuth(d::HorizontalDipole) = d.azimuth_angle
elevation(d::HorizontalDipole) = π/2
