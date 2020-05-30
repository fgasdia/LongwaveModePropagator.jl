"""
    Waveguide

A waveguide propagation path with a background `bfield`, ionosphere `constituent`,
and `ground`.
"""
abstract type Waveguide end

# TODO: constituent to species

"""
    HomogeneousWaveguide{B,C,G} <: Waveguide

Defines a homogeneous segment of waveguide with a `bfield`, `constituent`,
`ground`, and `distance` of the start of the segment from the `Emitter`.
"""
@with_kw struct HomogeneousWaveguide{B,C,G} <: Waveguide
    bfield::B  # BField
    constituent::C  # `Constituent` or `Vector{Constituent}`
    ground::G  # Ground
    distance::Float64
end
function HomogeneousWaveguide(bfield, constituent, ground) = HomogeneousWaveguide(bfield, constituent, ground, 0)

struct SegmentedWaveguide <: Waveguide
    v::Vector{HomogeneousWaveguide{F,G}}
end

# TODO: WKBWaveguide (that's why we can't call SegmentedWaveguide -> InhomogeneousWaveguide)
