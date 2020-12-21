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

# ## Building the mesh grid
#
# The U.S. Navy's Long Wavelength Propagation Capability searches the region of the
# complex plane from 30° to 90° in the real axis and 0° to -10° in the imaginary axis.
# The lowest order modes are nearest ``90° - i0°`` excluding 90°.
# Modes with large negative imaginary components are highly attenuated and have
# relatively little affect on the total field in the waveguide.
# Large negative imaginary angles also pose numerical issues for the computation of
# the modified Hankel functions of order one-third used to describe the height gain
# functions of the fields in the waveguide.

function defaultmesh(z_bl, z_tr, Δr)
    
end
