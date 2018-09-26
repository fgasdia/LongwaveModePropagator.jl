module PropagationPath

using StaticArrays

struct SegmentParameters
	dist::AbstractArray
	lat::AbstractArray
	lon::AbstractArray
	azim::AbstractArray
	dip::AbstractArray
	bfield::AbstractArray
	sigma::AbstractArray
	epsr::AbstractArray
	hprime::AbstractArray
	beta::AbstractArray
	numsegments::Integer
end

"""
    segmentpath(σ, ϵᵣ, h′, β, inputs)

Break propagation path into segments with appropriate magnetic field, ground, and ionosphere
parameters.

The LWPC LWP_PATH subroutine is called in place of preseg to automatically determine when a
new mode segment should be created. By default it traverses every 20km along the path and
has an algorithm that looks at the magnetic field, ground conductivity, and ionosphere
profile to determine if a new segment should be created or the previous segment extended.

This new function segmentpath() takes the place of LWP_PATH and LWP_PRESEG. It splits the
path up into equal intervals so we know how many intervals will be required a priori and
create the appropriately sized StaticArrays for each segment component. This will usually
mean more segments need to be computed than in typical LWPC and could result in longer run
times, but will result in a more "true" and smoother profile. The path components filled
include ionosphere profile, ground conductivity, and magnetic field.

Could pass dispatch differently on ground, bfield, and ionosphere. If ground is a float,
then assume it is constant everywhere. If ground is a function handle, then call that
function to get conductivity values (could be a wrapper to reading data file in a specific
format). If bfield is not supplied, then don't use bfield, if bfield is a 3 element vector, 
then assume constant bfield everywhere, if bfield is a function handle, then use it as a 
wrapper to reading datafile. If profile is a float, assume it is constant everywhere; it may
also be an array of [dist, h′, β], or a function handle.

Here is a `variable` in the code which represents LaTeX equation ``y = x + \\sin(2θ)``.

When necessary, provide an argument list.

# Arguments
- `n::Integer`: the number of elements to compute.
- `dim::Integer=1`: the dimensions along which to perform the computation.
"""
function segmentpath(σ::Float64, ϵᵣ::Float64, h′::Float64, β::Float64,
					 inputs::Inputs)
	numsegments = 1
	segments = SegmentParameters([0.0],
								 [inputs.transmitter.lat],
								 [inputs.transmitter.lon],
								 [missing],
								 [missing],
								 [missing],
								 [σ],
								 [ϵᵣ],
								 [h′],
								 [β],
								 numsegments)

	return segments
end

function segmentpath(σ::Float64, ϵᵣ::Float64,
					 dist::AbstractArray, h′::AbstractArray, β::AbstractArray,
					 inputs::Inputs)
	@assert size(h′) == size(β) == size(dist) "h′, β, and `dist` must be same size"

	# Use h′, β `dist` array as path segments
	numsegments = length(dist)
	segments = SegmentParameters(dist,
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 fill(missing, numsegments),
								 fill(missing, numsegments),
								 fill(missing, numsegments),
								 fill(σ, numsegments),
								 fill(ϵᵣ, numsegments),
								 h′,
								 β,
								 numsegments)

	for ii = 1:numsegments
		segments.lat[ii] = # Calculate
		segments.lon[ii] = # Calculate
	end

	return segments
end

function segmentpath(gnd::Function)

end

"""
The LWPC LWP_DRIVER subroutine "processes" the presegmented path parameters. Effectively
this means it loops over each segment, calls MF_DRIVER to get starting solutions, then SW/MF
as needed to find modes.

TODO: Figure out if the switching between MF and SW is necessary.
"""
function lwp_driver()
	numsegments = floor(Int, inputs.maxrange/inputs.deltarange)

	segments = SegmentParameters(Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 Array{Float64}(undef, numsegments),
								 numsegments)

	# First segment is the transmitter
	segments.dist[1] = 0.0
	segments.lat[1] = inputs.transmitter.lat
	segments.lon[1] = inputs.transmitter.lon
	segments.azim[1] = missing
	segments.dip[1] = missing
	segments.bfield[1] = missing
	segments.sigma[1] = σ
	segments.epsr[1] = ϵᵣ
	segments.hprime[1] = h′
	segments.beta[1] = β

	for ii = 2:numsegments

	end

	return segments
end

end