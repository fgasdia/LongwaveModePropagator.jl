# # Density and collision frequency as interpolating functions
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
# a linear interpolator and a cubic spline interpolator.
# We'll use [NormalHermiteSplines.jl](https://github.com/IgorKohan/NormalHermiteSplines.jl)
# to construct Hermite splines.
# 
# Unfortunately, Julia doesn't currently have a well-supported library with the
# [PCHIP](https://blogs.mathworks.com/cleve/2012/07/16/splines-and-pchips/) interpolator,
# so I must load `pchip` from my _private_ repository Interp.jl.

using Printf
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using Plots, Distances

using Interpolations, NormalHermiteSplines
using Interp

# Here's the discrete profile data. Some interpolators will benefit from random sample
# points, but real density data will almost always be on a grid.

zs = 0:500:110e3
Ne = waitprofile.(zs, 75, 0.32; cutoff_low=40e3)

# Let's construct the interpolators.

linear_itp = LinearInterpolation(zs, Ne)
cubic_itp = CubicSplineInterpolation(zs, Ne)

spline = prepare(collect(zs), RK_H1())
spline = construct(spline, Ne)
hermite_itp(z) = evaluate(spline, collect(z))

pchip_itp = pchip(zs, Ne)


# (a spline built with RK_H0 kernel is a continuous function,
#  a spline built with RK_H1 kernel is a continuously differentiable function,
#  a spline built with RK_H2 kernel is a twice continuously differentiable function).

zs_fine = 0:5:110e3
Ne_fine = waitprofile.(zs_fine, 75, 0.32; cutoff_low=40e3)

linear_fine = linear_itp.(zs_fine)
cubic_fine = cubic_itp.(zs_fine)
hermite_fine = hermite_itp(zs_fine)
pchip_fine = pchip_itp.(zs_fine)

# The profiles are compared using percentage difference relative to the true profile.

cmp(a,b) = (a .- b)./b.*100

dNe = cmp(Ne_fine, Ne_fine)
dlinear = cmp(linear_fine, Ne_fine)
dcubic = cmp(cubic_fine, Ne_fine)
dhermite = cmp(hermite_fine, Ne_fine)
dpchip = cmp(pchip_fine, Ne_fine)

# To plot densities with a log scale, we set values less than 0.1 to NaN.

cl(x) = replace(v->v <= 0.1 ? NaN : v, x)
lc(x) = replace(x, NaN => 0)

p1 = plot(cl([Ne_fine linear_fine cubic_fine hermite_fine pchip_fine]),
    zs_fine/1000, xscale=:log10, xlabel="Ne (m⁻³)", ylabel="Altitude (km)",
    legend=:topleft, labels=["Truth" "Linear" "Cubic" "Hermite" "PCHIP"])
p2 = plot(lc([dNe dlinear dcubic dhermite dpchip]),
    zs_fine/1000, xlabel="% difference", legend=false, xlims=(-100, 100))
p3 = plot(lc([dNe dlinear dcubic dhermite dpchip]),
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

interpolators = 

for (n, v) in ("linear"=>linear_fine, "cubic"=>cubic_fine, "hermite"=>hermite_fine, "pchip"=>pchip_fine)
    @printf("%s: %.3e\n",n, rmsd(v, Ne_fine))
end

# ## Propagation results
# 
# Let's compare propagation results and performance with each of these interpolators.


# Here are the `Species` for each interpolator.
# They will all use the analytic [`electroncollisionfrequency`](@ref).

interpolators = (
    "truth" => z->waitprofile(z, 75, 0.32; cutoff_low=40e3),
    "linear"=> linear_itp,
    "cubic" => cubic_itp,
    "hermite" => hermite_itp,
    "pchip" => pchip_itp
)

function propagateitp(interpolators)
    bfield = BField(50e-6, deg2rad(68), deg2rad(111))
    ground = Ground(5, 0.00005)

    tx = Transmitter(24e3)
    rx = GroundSampler(0:5e3:3000e3, Fields.Ez)

    results = Dict{String,NTuple{3,Float64}}()
    for (n, itp) in interpolators
        species = Species(QE, ME, itp, electroncollisionfrequency)
        waveguide = HomogeneousWaveguide(bfield, electrons, ground)
        E, amp, phase = propagate(waveguide, tx, rx)
    end
    return results
end


# Timing can vary when run by github actions build, but normalized results from a local run are here:

p1 = plot(rx.distance/1000, ae, label="electrons", ylabel="Amplitude (dB μV/m)")
plot!(p1, rx.distance/1000, aei, label="electrons & ions")
p2 = plot(rx.distance/1000, aei-ae,
    ylims=(-0.5, 0.5), xlabel="Range (km)", ylabel="Δ", legend=false)
plot(p1, p2, layout=grid(2,1,heights=[0.7, 0.3]))
