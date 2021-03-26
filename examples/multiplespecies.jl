# # Multiple ionospheric species
# 
# This example demonstrates an ionosphere with multiple constituents -- not only electrons,
# but ions as well.
# 
# ## Background
# 
# As discussed in the background of [Interpreting h′ and β](@ref interpreting-hpbeta),
# the propagation of VLF waves depends on the number density and collision frequency of
# charged species. Together, these two quantities describe the conductivity profile of the
# ionosphere. Not only electrons, but massive ions, influence the conductivity of the
# D-region. The effect of multiple species are combined simply by summing the susceptibility
# tensor elements for each species' density and collision frequency. 
# 
# ## Implementation
# 
# There's nothing special about electrons in LongwaveModePropagator.
# In all of the other examples, the only ionospheric constituent is electrons, but they
# are represented as a [`Species`](@ref) type.
# LongwaveModePropagator let's you define multiple ionospheric `Species` and pass them as a
# `tuple` to most functions that take the argument `species`.
# 
# !!! note
#   Multiple species should be passed as a `tuple` of `Species` for performance.
#   Although multiple species can also be passed as a `Vector`, the type will be
#   `Vector{Species}`. This is not a concrete type (because `Species` are parameterized on
#   their `numberdensity` and `collisionfrequency` functions) and the compiler will not emit
#   performant code.


# Let's compare electrons-only and multiple-species ionospheres.

using Plots
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME

# `QE` is literally the charge of an electron and is therefore negative.

QE
