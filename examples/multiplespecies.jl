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

# The [`Species`](@ref) type defines the `charge`, `mass`, the `numberdensity` as a function
# of altitude, and the neutral `collisionfrequency` as a function of altitude.
# 
# We'll import the internally defined electron charge `QE` and mass `ME` for convenience.

using LongwaveModePropagator: QE, ME

# `QE` is literally the charge of an electron and is therefore negative.

QE

# LMP also exports the function [`waitprofile`](@ref) for the Wait and Spies (1964) electron
# density, [`electroncollisionfrequency`](@ref) for the electron-neutral collision frequency,
# and [`ioncollisionfrequency`](@ref) for the ion-neutral collision frequency.
# 
# LWPC supports built-in electrons, a positive ion, and a negative ion.
# Each of the ions have a mass of 58,000 electrons (approximately the mass of O₂).
# 
# Here we'll use the `waitprofile` for the electrons, but make up a not-very-realistic
# exponential profile for the positive ion.

Ne(z) = waitprofile(z, 75, 0.32)
Np(z) = 2e6*exp(1e-4*z)

# We define density of the negative ions to conserve charge neutrality. 

Nn(z) = Np(z) - Ne(z)

# For plotting, we'll replace densities below 1 with `NaN`.

mask(x) = x < 1 ? NaN : x

z = 0:1e3:110e3
plot(mask.(Ne.(z)), z/1000, label="Ne", linewidth=1.5)
plot!(mask.(Np.(z)), z/1000, label="Np", linewidth=1.5)
plot!(mask.(Nn.(z)), z/1000, label="Nn", linewidth=1.5)
plot!(xlabel="Density (m⁻³)", ylabel="Altitude (km)",
    xscale=:log10, legend=:topleft)
#md savefig("multiplespecies_density.png"); nothing # hide
#md # ![](multiplespecies_density.png)

# Here are the electron and ion collision frequencies:

plot(electroncollisionfrequency.(z), z/1000, label="νe", linewidth=1.5)
plot!(ioncollisionfrequency.(z), z/1000, label="νi", linewidth=1.5)
plot!(xlabel="Collision frequency (s⁻¹)", ylabel="Altitude (km)", xscale=:log10)
#md savefig("multiplespecies_collisionfreq.png"); nothing # hide
#md # ![](multiplespecies_collisionfreq.png)

# The [`Species`](@ref):

electrons = Species(QE, ME, Ne, electroncollisionfrequency)
posions = Species(abs(QE), 58000*ME, Np, ioncollisionfrequency)
negions = Species(QE, 58000*ME, Nn, ioncollisionfrequency)

# We'll compare an electrons-only and electrons-ions ionosphere.

tx = Transmitter(24e3)
rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

bfield = BField(50e-6, π/2, 0)
ground = GROUND[10]

ewvg = HomogeneousWaveguide(bfield, electrons, ground)
eiwvg = HomogeneousWaveguide(bfield, (electrons, posions, negions), ground)

Ee, ae, pe = propagate(ewvg, tx, rx)
Eei, aei, pei  = propagate(eiwvg, tx, rx)

p1 = plot(rx.distance/1000, ae, label="electrons", ylabel="Amplitude (dB μV/m)")
plot!(p1, rx.distance/1000, aei, label="electrons & ions")
p2 = plot(rx.distance/1000, aei-ae, xlabel="Range (km)", ylabel="Δ", legend=false)
plot(p1, p2, layout=grid(2,1,heights=[0.7, 0.3]))
#md savefig("multiplespecies_amplitude.png"); nothing # hide
#md # ![](multiplespecies_amplitude.png)


# NOTE ABOUT TIMING (BOTH BETWEEN E AND EI AND VECTOR VS TUPLE)