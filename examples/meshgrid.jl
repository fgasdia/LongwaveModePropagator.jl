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
pkg"add LongwaveModePropagator, Plots"

# Then we import the packages used in this example.

using Plots
using LongwaveModePropagator
using LongwaveModePropagator: solvemodalequation

## Using the GR backend
gr(legend=false)

# ## Building the mesh grid
#
# The U.S. Navy's Long Wavelength Propagation Capability searches the region of the
# complex plane from 30° to 90° in the real axis and 0° to -10° in the imaginary axis.
# The lowest order modes are nearest ``90° - i0°``, excluding 90°.
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

function rectanglemesh(zbl, ztr, Δr)
    rzbl, izbl = reim(zbl)
    rztr, iztr = reim(ztr)

    X = rztr - rzbl
    Y = iztr - izbl

    n = ceil(Int, Y/Δr)
    dy = Y/n

    # dx = sqrt(Δr² - (dy/2)²), solved for equilateral triangle
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m
    half_dx = dx/2

    T = promote_type(ComplexF64, typeof(zbl), typeof(ztr), typeof(Δr))
    v = Vector{T}()

    shift = false  # we will displace every other line
    for j = 0:n
        y = izbl + dy*j

        for i = 0:m
            x = rzbl + dx*i

            if shift
                x -= half_dx
            end

            if (i ) == m
                shift = !shift
            end

            push!(v, complex(x, y))
        end
    end

    return v
end

zbl = complex(30.0, -10.0)
ztr = complex(89.9, 0.0)
Δr = 0.5
v = rectanglemesh(zbl, ztr, Δr)

# Let's plot `v` to make sure it's correct.
# We set the axis limits so there are equal ranges in both axes.
# The maximum real and imaginary values are drawn in red

plot(real(v), imag(v), seriestype=:scatter,
     xlims=(80, 90), ylims=(-10, 0),
     xlabel="real(θ)", ylabel="imag(θ)")
plot!([80, 90], [0, 0], color="red")
plot!([89.9, 89.9], [-10, 0], color="red")

# Looks good!
#
# Naturally we want to be able to use the coarsest initial mesh grid as possible
# because the reflection coefficient has to be integrated through the ionosphere
# for every node of the grid.
# But for now, let's solve the modal equation on a very fine grid to explore its
# phase across the complex quadrant.


function defaultmesh(z_bl, z_tr, Δr)

end
