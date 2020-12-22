# # Mesh grid for mode finding
#
# The [Global complex Roots and Poles Finding](https://doi.org/10.1109/TAP.2018.2869213)
# (GRPF) algorithm searches for roots and poles of complex-valued functions by sampling
# the function at the nodes of a regular mesh and analyzing the complex phase.
# The mesh search region can be arbitrarily shaped, and guarantees that none of
# the roots/poles can be missed as long as the mesh step is sufficiently small.
# Specifically, the phase change of the function value cannot exceed three quadrants
# across mesh edges. Because the GRPF algorithm uses an adaptive meshing process,
# the initial mesh node spacing can actually be greater than the distance between
# the roots/poles of the function, but there is no established method for _a priori_
# estimation of the required initial sampling of an arbitrary function.
# The GRPF algorithm is used by `LongwaveModePropagator.jl` to identify eigenangles
# of `HomogeneousWaveguide`'s.
# Therefore, we experimentally determine a sufficient initial mesh grid spacing to
# identify roots of the mode equation.
#
# ## Mode equation
#
# The criteria for resonance in the Earth-ionosphere waveguide is described by the
# fundamental equation of mode theory:
#
# ```math
# \det(\overline{\bm{R}}(\theta)\bm{R}(\theta) - \bm{I}) = 0
# ```
#
# where ``\overline{\bm{R}}(\theta)`` is the reflection coefficient matrix for the
# ground and ``\bm{R}(\theta)`` is the reflection coefficient of the ionosphere.
# Both are functions of the complex angle ``\theta``. The discrete values of
# ``\theta`` for which the fundamental equation is satisfied are known as
# _eigenangles_.
#
# ## Install dependencies
#
# First we make sure all required packages are installed.

using Pkg
## Pkg.add("LongwaveModePropagator")
pkg"add Plots, DisplayAs, RootsAndPoles, VoronoiDelaunay"

# Then we import the packages used in this example.

using Plots, DisplayAs
using RootsAndPoles

using ..LongwaveModePropagator
using ..LongwaveModePropagator: solvemodalequation, QE, ME

## Using the GR backend
gr(legend=false);

# ## Building the mesh grid
#
# The U.S. Navy's Long Wavelength Propagation Capability searches the region of the
# complex plane from ``30°`` to ``90°`` in the real axis and ``0°`` to ``-10°`` in the
# imaginary axis.
# The lowest order modes are nearest ``90° - i0°``, excluding ``90°``.
# Modes with large negative imaginary components are highly attenuated and have
# relatively little affect on the total field in the waveguide.
# Large negative imaginary angles also pose numerical issues for the computation of
# the modified Hankel functions of order one-third used to describe the height gain
# functions of the fields in the waveguide.
#
# GRPF is most efficient when the initial mesh grid consists of equilateral triangles.
# We can construct such a grid with the function below, which builds the grid over a
# rectangular region from `zbl` in the bottom left to `ztr` in the top right
# of the lower right quadrant of the complex plane (with positive real and negative
# imaginary) and `Δr` is the mesh step.
#

function rectanglemesh(zbl, ztr, Δr)
    rzbl, izbl = reim(zbl)
    rztr, iztr = reim(ztr)

    X = rztr - rzbl
    Y = iztr - izbl

    n = ceil(Int, Y/Δr)
    dy = Y/n

    ## dx = sqrt(Δr² - (dy/2)²), solved for equilateral triangle
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m
    half_dx = dx/2

    T = promote_type(ComplexF64, typeof(zbl), typeof(ztr), typeof(Δr))
    v = Vector{T}()

    shift = false  # we will displace every other line by dx/2
    for j = 0:n
        y = izbl + dy*j

        for i = 0:m
            x = rzbl + dx*i

            if shift
                x -= half_dx
            end

            if i == m
                shift = !shift
            end

            push!(v, complex(x, y))
        end
    end

    return v, (dy, dx)
end

zbl = complex(30.0, -10.0)
ztr = complex(89.9, 0.0)
Δr = 0.5
v, dydx = rectanglemesh(zbl, ztr, Δr);

# Let's plot `v` to make sure it's correct.
# We set the axis limits so there are equal ranges in both axes.
# The maximum real and imaginary values are drawn with red lines.

img = plot(real(v), imag(v), seriestype=:scatter,
           xlims=(80, 90), ylims=(-10, 0),
           xlabel="real(θ)", ylabel="imag(θ)",
           size=(450,375))
plot!(img, [80, 90], [0, 0], color="red")
plot!(img, [89.9, 89.9], [-10, 0], color="red")
DisplayAs.PNG(img)  #hide

# Looks good!
#
# ## Exploring the modal equation phase
#
# Naturally we want to be able to use the coarsest initial mesh grid as possible
# because the reflection coefficient has to be integrated through the ionosphere
# for every node of the grid.
# But for now, let's solve the modal equation on a fine grid to explore its
# phase across a larger region of the bottom right complex quadrant.
#
# We'll define four different `PhysicalModeEquation` objects from
# `LongwaveModePropagator.jl` that we'll use throughout this example.

lowfrequency = Frequency(10e3)
midfrequency = Frequency(20e3)
highfrequency = Frequency(100e3)

day = Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)
night = Species(QE, ME, z->waitprofile(z, 85, 0.9), electroncollisionfrequency)

daywaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), day, GROUND[5])
nightwaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), night, GROUND[5])

day_low_me = PhysicalModeEquation(lowfrequency, daywaveguide)
day_mid_me = PhysicalModeEquation(midfrequency, daywaveguide)
day_high_me = PhysicalModeEquation(highfrequency, daywaveguide)

night_mid_me = PhysicalModeEquation(midfrequency, nightwaveguide);

day_low_title = "10 kHz\nh′: 75, β: 0.35"
day_mid_title = "20 kHz\nh′: 75, β: 0.35"
day_high_title = "100 kHz\nh′: 75, β: 0.35"
night_mid_title = "20 kHz\nh′: 85, β: 0.9";

# Then construct the mesh.

zbl = complex(0.0, -40.0)
ztr = complex(89.9, 0.0)
Δr = 0.2
v, dydx = rectanglemesh(zbl, ztr, Δr);

# Now we simply iterate over each node of the mesh, evaluating the modal equation with
# `solvemodalequation` from `LongwaveModePropagator.jl`.
# We use Julia's `Threads.@threads` multithreading capability to speed up the computation.

function modeequationphase(me, mesh)
    phase = Vector{Float64}(undef, length(mesh))
    Threads.@threads for i in eachindex(mesh)
        f = solvemodalequation(deg2rad(mesh[i]), me)
        phase[i] = rad2deg(angle(f))
    end
    return phase
end

phase = modeequationphase(day_mid_me, v);

# Plotting the results, we can see that there are clearly identifiable locations where
# white, black, blue, and orange all meet.
# Each of these locations are either a root or pole in the daytime ionosphere.

yvec = imag(zbl):dydx[1]:imag(ztr)
xvec = real(zbl):dydx[2]:real(ztr)
img = heatmap(xvec, yvec, reshape(phase, length(xvec), length(yvec))',
              color=:twilight, clims=(-180, 180),
              xlims=(0, 90), ylims=(-40, 0),
              xlabel="real(θ)", ylabel="imag(θ)",
              legend=true, size=(600, 400),
              title=day_mid_title)
DisplayAs.PNG(img)  #hide

# If we switch to a nighttime ionosphere with a high Wait β parameter, we see that the
# roots/poles move closer to the axes.
# A perfectly reflecting conductor has eigenangles along the real and complex axes.

phase = modeequationphase(night_mid_me, v);

img = heatmap(xvec, yvec, reshape(phase, length(xvec), length(yvec))',
              color=:twilight, clims=(-180, 180),
              xlims=(0, 90), ylims=(-40, 0),
              xlabel="real(θ)", ylabel="imag(θ)",
              legend=true, size=(600, 400),
              title=night_mid_title)
DisplayAs.PNG(img)  #hide

# At lower frequencies, the roots/poles move further apart.

phase = modeequationphase(day_low_me, v);

img = heatmap(xvec, yvec, reshape(phase, length(xvec), length(yvec))',
              color=:twilight, clims=(-180, 180),
              xlims=(0, 90), ylims=(-40, 0),
              xlabel="real(θ)", ylabel="imag(θ)",
              legend=true, size=(600, 400),
              title=day_low_title)
DisplayAs.PNG(img)  #hide

# ## Modifying the mesh grid
#
# No roots or poles appear in the lower right triangle of the domain.
# Even if they did, they would be highly attenuated modes.
# Therefore, to save compute time, we can exclude the lower right triangle of the domain
# from the initial mesh.
# The function is simply the `rectanglemesh` with an extra check to see if the
# node is in the upper left triangle or not.
# One important difference is that units now matter.
# Although `rectanglemesh` works where `zbl`, `ztr`, and `Δr` are in degrees or radians,
# `π/2` is hardcoded into `trianglemesh` and therefore the arguments must also be in
# radians.

function trianglemesh(zbl, ztr, Δr)
    rzbl, izbl = reim(zbl)
    rztr, iztr = reim(ztr)

    X = rztr - rzbl
    Y = iztr - izbl

    n = ceil(Int, Y/Δr)
    dy = Y/n

    ## dx = sqrt(Δr² - (dy/2)²), solved for equilateral triangle
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m
    half_dx = dx/2

    slope = 1  # 45° angle

    T = promote_type(ComplexF64, typeof(zbl), typeof(ztr), typeof(Δr))
    v = Vector{T}()

    shift = false  # we will displace every other line by dx/2
    for j = 0:n
        y = izbl + dy*j

        for i = 0:m
            x = rzbl + dx*i

            if shift
                x -= half_dx
            end

            if i == m
                shift = !shift
            end

            ## NEW: check if `x, y` is in upper left triangle
            if y >= slope*x - π/2
                push!(v, complex(x, y))
            end
        end
    end

    return v, (dy, dx)
end

zbl = deg2rad(complex(30.0, -10.0))
ztr = deg2rad(complex(89.9, 0.0))
Δr = deg2rad(0.5)

v, dydx = trianglemesh(zbl, ztr, Δr);

# We convert back to degrees just for plotting.

vdeg = rad2deg.(v)

img = plot(real(vdeg), imag(vdeg), seriestype=:scatter,
           xlims=(80, 90), ylims=(-10, 0),
           xlabel="real(θ)", ylabel="imag(θ)",
           size=(450,375))
plot!(img, [80, 90], [0, 0], color="red")
plot!(img, [0, 90], [-90, 0], color="red")
DisplayAs.PNG(img)  #hide

# The example continues in [Mesh grid for mode finding - Part 2](@ref).
