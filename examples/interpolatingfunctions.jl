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

using Interpolations, NormalHermiteSplines
using Interp

using Plots

# Here's the discrete profile data. Some interpolators will benefit from random sample
# points, but real density data will almost always be on a grid.

zs = 0:500:110e3
Ne = waitprofile.(zs, 75, 0.32; cutoff_low=40e3)

# Let's construct the interpolators.

linear_itp = LinearInterpolation(zs, Ne)
cubic_itp = CubicSplineInterpolation(zs, Ne)

spline = prepare(collect(zs), RK_H1())
spline = construct(spline, Ne)
hermite_itp(z) = evaluate(spline, z)

pchip_itp = pchip(zs, Ne)


# (a spline built with RK_H0 kernel is a continuous function,
    #  a spline built with RK_H1 kernel is a continuously differentiable function,
    #  a spline built with RK_H2 kernel is a twice continuously differentiable function).

zs_fine = 0:5:110e3
Ne_fine = waitprofile.(zs_fine, 75, 0.32; cutoff_low=40e3)

# To plot correctly with log10 xscale, set NaN to NaN

cl(x) = replace(v->v <= 0.1 ? NaN : v), x)

plot(cl([Ne_fine linear_itp.(zs_fine) cubic_itp.(zs_fine) hermite_itp.(zs_fine) pchip_itp.(zs_fine)]),
    zs_fine/1000, xscale=:log10, labels=["Truth" "Linear" "Cubic" "Hermite" "PCHIP"])


# ## Propagation results
# 
# and performance... (table?)

# Here are the `Species` for each interpolator.

linear_species = Species(QE, ME, linear_itp, electroncollisionfrequency)
cubic_species = Species(QE, ME, cubic_itp, electroncollisionfrequency)
hermite_species = Species(QE, ME, hermite, electroncollisionfrequency)
pchip_species = Species(QE, ME, pchip_itp, electroncollisionfrequency)

bfield = BField(50e-6, deg2rad(68), deg2rad(111))
ground=Ground(5, 0.00005)

tx = Transmitter(24e3)
rx = GroundSampler(0:5e3:3000e3, Fields.Ez)

waveguide = HomogeneousWaveguide(bfield, species, ground)

E, amp, phase = propagate(waveguide, tx, rx)
