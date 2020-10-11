"""
    Waveguide

A waveguide propagation path with a background `BField`, ionosphere `Species`,
and `Ground`.
"""
abstract type Waveguide end

"""
    HomogeneousWaveguide{S} <: Waveguide

Defines a homogeneous segment of waveguide with a `bfield::BField`, `species::S`,
`ground::Ground`, and `distance` of the start of the segment from the `Emitter`.
"""
struct HomogeneousWaveguide{S} <: Waveguide
    bfield::BField
    species::S  # `Species` or `Vector{Species}`
    ground::Ground
    distance::Float64
end
HomogeneousWaveguide(bfield, species::S, ground) where S =
    HomogeneousWaveguide{S}(bfield, species, ground, 0.0)

struct SegmentedWaveguide{T<:Vector{<:Waveguide}} <: Waveguide
    v::T
end
Base.length(w::SegmentedWaveguide) = length(w.v)
Base.sort!(w::SegmentedWaveguide) = sort!(w.v, by=x->getfield(x, :distance))
Base.sort(w::SegmentedWaveguide) = SegmentedWaveguide(sort!(copy(w.v)))
Base.issorted(w::SegmentedWaveguide) = issorted(w.v, by=x->getfield(x, :distance))
Base.size(w::SegmentedWaveguide) = size(w.v)
Base.getindex(w::SegmentedWaveguide, i::Int) = w.v[i]
Base.setindex(w::SegmentedWaveguide, x, i::Int) = (w.v[i] = x)
Base.push!(w::SegmentedWaveguide, x) = push!(w.v, x)

"""
    SegmentedWaveguide(w<:Waveguide)

Convenience function that generates a `SegmentedWaveguide` with an empty `Vector`
of type `w`.
"""
SegmentedWaveguide(::Type{T}) where T<:Waveguide = SegmentedWaveguide(Vector{T}())

# TODO: WKBWaveguide (that's why we can't call SegmentedWaveguide -> InhomogeneousWaveguide)
