"""
    Waveguide

A waveguide propagation path with a background `BField`, ionosphere `Species`, and `Ground`.
"""
abstract type Waveguide end

"""
    HomogeneousWaveguide{S} <: Waveguide

Defines a homogeneous segment of waveguide.

# Fields

- `bfield::BField`: background magnetic field.
- `species::S`: ionosphere constituents.
- `ground::Ground`: waveguide ground.
- `distance::Float64`: distance from the `Emitter` at the start of the segment in meters.
"""
struct HomogeneousWaveguide{S} <: Waveguide
    bfield::BField
    species::S  # `Species` or `Vector{Species}`
    ground::Ground
    distance::Float64
end

"""
    HomogeneousWaveguide(bfield, species, ground)

By default, `distance` is `0.0`.
"""
HomogeneousWaveguide(bfield, species::S, ground) where S =
    HomogeneousWaveguide{S}(bfield, species, ground, 0.0)

"""
    adjoint(w::HomogeneousWaveguide)

Return `w` with an adjoint `BField` having an `x` component of opposite sign.
"""
function adjoint(w::HomogeneousWaveguide)
    @unpack bfield, species, ground, distance = w
    adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
    return HomogeneousWaveguide(adjoint_bfield, species, ground, distance)
end

"""
    SegmentedWaveguide{T<:Vector{<:Waveguide}} <: Waveguide

A collection of `Waveguide`s make up an inhomogeneous segmented waveguide.
"""
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
Base.iterate(w::SegmentedWaveguide) = iterate(w.v)
Base.iterate(w::SegmentedWaveguide, state) = iterate(w.v, state)

# TODO: WKBWaveguide (that's why we can't call SegmentedWaveguide -> InhomogeneousWaveguide)
