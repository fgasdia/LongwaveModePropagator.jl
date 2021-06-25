# # [Density and collision frequency as interpolating functions](@id interpolating-functions)
# 
# The ionospheric [`Species`](@ref) type used when defining a [`HomogeneousWaveguide`](@ref)
# requires `numberdensity` and `collisionfrequency` fields to be callable functions of
# altitude `z` in meters that return `Float64` values of number density in m⁻³ and 
# species-neutral collision frequency in s⁻¹.
# 
# Sometimes a simple analytical function can be provided, such as [`waitprofile`](@ref).
# Other times, the profile is constructed from tabular data at discrete heights which must
# be wrapped in an interpolator before being passed to `Species`.
# 
# The choice of interpolating function not only impacts the accuracy of the propagation
# result (how well the interpolation resembles the true profile), but it can also impact
# the runtime of [`propagate`](@ref) and related functions.
# 
# In this example we'll compare a few different interpolating functions.  

# ## Profiles
# 
# We'll use a [`waitprofile`](@ref) with a sharp cutoff at 40 km to compare the
# interpolating functions, but more complex profiles may not be strictly monotonically
# increasing.
# 
# From [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) we'll use
# a linear interpolator, a cubic spline interpolator, and monotonic interpolators.
# We'll use [NormalHermiteSplines.jl](https://github.com/IgorKohan/NormalHermiteSplines.jl)
# to construct Hermite splines.

using Printf
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using Plots, Distances

using Interpolations, NormalHermiteSplines
nothing #hide

# Here's the discrete profile data. Some interpolators will benefit from random sample
# points, but most density data will be on a grid.

zs = 0:1e3:110e3
Ne = waitprofile.(zs, 75, 0.32; cutoff_low=40e3);

# Let's construct the interpolators.

linear_itp = LinearInterpolation(zs, Ne)
cubic_itp = CubicSplineInterpolation(zs, Ne);

# There's not great documentation on the monotonic interpolators of Interpolations.jl as of
# `v0.13`, but several are supported.

fb_itp = interpolate(zs, Ne, FritschButlandMonotonicInterpolation())
fc_itp = interpolate(zs, Ne, FritschCarlsonMonotonicInterpolation())
s_itp = interpolate(zs, Ne, SteffenMonotonicInterpolation());

# And the Hermite splines

spline = prepare(collect(zs), RK_H1())
spline = construct(spline, Ne)
hermite_itp(z) = evaluate_one(spline, z);

# (a spline built with RK_H0 kernel is a continuous function,
#  a spline built with RK_H1 kernel is a continuously differentiable function,
#  a spline built with RK_H2 kernel is a twice continuously differentiable function).

zs_fine = 0:5:110e3
Ne_fine = waitprofile.(zs_fine, 75, 0.32; cutoff_low=40e3)

linear_fine = linear_itp.(zs_fine)
cubic_fine = cubic_itp.(zs_fine)
fb_fine = fb_itp.(zs_fine)
fc_fine = fc_itp.(zs_fine)
s_fine = s_itp.(zs_fine)
hermite_fine = hermite_itp.(zs_fine);

# The profiles are compared using percentage difference relative to the true profile.

cmp(a,b) = (a .- b)./b.*100

dNe = cmp(Ne_fine, Ne_fine)
dlinear = cmp(linear_fine, Ne_fine)
dcubic = cmp(cubic_fine, Ne_fine)
dfb = cmp(fb_fine, Ne_fine)
dfc = cmp(fc_fine, Ne_fine)
ds = cmp(s_fine, Ne_fine)
dhermite = cmp(hermite_fine, Ne_fine);

# To plot densities with a log scale, we set values less than 0.1 to NaN.

cl(x) = replace(v->v <= 0.1 ? NaN : v, x)
lc(x) = replace(x, NaN => 0)

p1 = plot(cl([Ne_fine linear_fine cubic_fine hermite_fine fb_fine fc_fine s_fine]),
    zs_fine/1000, xscale=:log10, xlabel="Ne (m⁻³)", ylabel="Altitude (km)",
    legend=:topleft, labels=["Truth" "Linear" "Cubic" "Hermite" "FritschButland" "FritschCarlson" "Steffen"])
p2 = plot(lc([dNe dlinear dcubic dhermite dfb dfc ds]),
    zs_fine/1000, xlabel="% difference", legend=false, xlims=(-100, 100))
p3 = plot(lc([dNe dlinear dcubic dhermite dfb dfc ds]),
    zs_fine/1000, xlabel="% difference", legend=false, xlims=(-0.1, 0.1))
plot(p1, p2, p3, layout=(1,3), size=(800,400))
#md savefig("interpolatingfunctions_profiles.png"); nothing # hide
#md # ![](interpolatingfunctions_profiles.png)
    
# Unsurprisingly, the error is highest at the cutoff altitude of 40 km where the densities
# below are 0.
# The linear interpolation has positive-biased errors at all other heights because linear
# interpolation does not accurately capture the true exponential profile.
# To avoid this, the interpolator could have been constructed with a finer initial grid, but
# that is not always possible.
# 
# The total RMS differences between each interpolator and the truth are:

for (n, v) in ("linear"=>linear_fine, "cubic"=>cubic_fine, "hermite"=>hermite_fine,
    "FritschButland"=>fb_fine, "FritschCarlson"=>fc_fine, "Steffen"=>s_fine)
    @printf("%s: %.3e\n",n, rmsd(v, Ne_fine))
end

# It's also important to note that the results may be different for non-exponential
# ionospheres with more complicated profiles. 
# 
# ## Propagation results
# 
# Let's compare propagation results and performance with each of these interpolators.
# Unfortunately, as of `v0.5.1` of NormalHermiteSplines.jl, the implementation of the
# package isn't efficient (see #3)[https://github.com/IgorKohan/NormalHermiteSplines.jl/issues/3]
# and because `hermite_itp` would be called millions of times in `propagate`, it is
# prohibitively slow.
# 
# Here are the `Species` for each of the other interpolators.
# They will all use the analytic [`electroncollisionfrequency`](@ref).

interpolators = (
    "truth" => z->waitprofile(z, 75, 0.32; cutoff_low=40e3),
    "linear"=> linear_itp,
    "cubic" => cubic_itp,
    "FritschButland" => fb_itp,
    "FritschCarlson" => fc_itp,
    "Steffen" => s_itp
)

function propagateitp(interpolators)
    bfield = BField(50e-6, deg2rad(68), deg2rad(111))
    ground = Ground(5, 0.00005)

    tx = Transmitter(24e3)
    rx = GroundSampler(0:5e3:3000e3, Fields.Ez)

    results = Dict{String,Tuple{Float64,Vector{Float64},Vector{Float64}}}()
    for (n, itp) in interpolators
        species = Species(QE, ME, itp, electroncollisionfrequency)
        waveguide = HomogeneousWaveguide(bfield, species, ground)
        t0 = time()
        _, amp, phase = propagate(waveguide, tx, rx)
        runtime = time() - t0
        results[n] = (runtime, amp, phase)
    end
    return results
end

propagateitp(interpolators);  # warmup
results = propagateitp(interpolators);

# First let's evaluate the propagation results.

d = 0:5:3000
p1 = plot(ylabel="Amplitude (dB μV/m)")
p2 = plot(xlabel="Range (km)", ylims=(-0.03, 0.03), ylabel="Δ", legend=false)
for (n, v) in results
    plot!(p1, d, v[2], label=n)
    plot!(p2, d, v[2]-results["truth"][2])
end
plot(p1, p2, layout=grid(2,1,heights=[0.7, 0.3]))
#md savefig("interpolatingfunctions_amp.png"); nothing # hide
#md # ![](interpolatingfunctions_amp.png)

p1 = plot(ylabel="Phase (deg)")
p2 = plot(xlabel="Range (km)", ylims=(-0.4, 0.4), ylabel="Δ", legend=false)
for (n, v) in results
    plot!(p1, d, rad2deg.(v[3]), label=n)
    plot!(p2, d, rad2deg.(v[3])-rad2deg.(results["truth"][3]))
end
plot(p1, p2, layout=grid(2,1,heights=[0.7, 0.3]))
#md savefig("interpolatingfunctions_phase.png"); nothing # hide
#md # ![](interpolatingfunctions_phase.png)

# Here is the mean absolute amplitude difference between each technique and the true profile:

for n in ("linear", "cubic", "FritschButland", "FritschCarlson", "Steffen")
    @printf("%s: %.3e\n",n, meanad(results[n][2], results["truth"][2]))
end

# The amplitude and phase for each of the interpolators matches the true exponential profile
# extremely closely.
#  
# Timing can vary when run by GitHub to build the documentation, so results here are from a
# local run:
# 
# |  Interpolator  | Runtime relative to `truth` |
# | -------------- | --------------------------- |
# | truth          |              1              |
# | linear         |             1.35            |
# | cubic          |             1.11            |
# | FritschButland |             1.23            |
# | FritschCarlson |             1.07            |
# | Steffen        |             1.26           |
