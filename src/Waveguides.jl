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

struct SegmentedWaveguide{T<:Waveguide} <: Waveguide
    v::Vector{T}
end
Base.length(w::SegmentedWaveguide) = length(w.v)
Base.sort!(w::SegmentedWaveguide) = sort!(w.v, by=x->getfield(x, :distance))
Base.sort(w::SegmentedWaveguide) = SegmentedWaveguide(sort!(copy(w.v)))
Base.issorted(w::SegmentedWaveguide) = issorted(w.v, by=x->getfield(x, :distance))
Base.getindex(w::SegmentedWaveguide, i::Int) = w.v[i]



# TODO: WKBWaveguide (that's why we can't call SegmentedWaveguide -> InhomogeneousWaveguide)
