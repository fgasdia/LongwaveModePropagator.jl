"""
    Waveguide

A waveguide propagation path with a background `BField`, ionosphere `Species`,
and `Ground`.
"""
abstract type Waveguide end

"""
    HomogeneousWaveguide{B,C,G} <: Waveguide

Defines a homogeneous segment of waveguide with a `BField`, `Species`,
`Ground`, and `distance` of the start of the segment from the `Emitter`.
"""
@with_kw struct HomogeneousWaveguide{B,S,G} <: Waveguide
    bfield::B  # BField
    species::S  # `Species` or `Vector{Species}`
    ground::G  # Ground
    distance::Float64
end
HomogeneousWaveguide(bfield, species, ground) = HomogeneousWaveguide(bfield, species, ground, 0)

struct SegmentedWaveguide <: Waveguide
    v::Vector{HomogeneousWaveguide}  # XXX: better typing? but not really necessary...
end

# TODO: WKBWaveguide (that's why we can't call SegmentedWaveguide -> InhomogeneousWaveguide)
