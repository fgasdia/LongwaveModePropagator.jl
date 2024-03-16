# # Ground
# 
# LongwaveModePropagator treats the ground boundary of the Earth-ionosphere waveguide
# as a homogeneous, nonmagnetic material described by its electrical conductivity ``\sigma``
# and relative permittivity ``\epsilon_r``.
# 
# In this example we'll compare the effect of ground conductivity on amplitude along
# the ground in the waveguide for day and nighttime ionospheres.

using Plots, Printf
using LongwaveModePropagator
const LMP = LongwaveModePropagator
nothing  #hide

# LWPC uses standard ground indices to describe combinations of relative permittivity
# and conductivity.
# The standard indices can be accessed from LMP as [`GROUND`](@ref).

sort(GROUND)

# We'll define all of the constant types and the day and night profiles.

## Constant values
const BFIELD = BField(50e-6, π/2, 0)

const TX = Transmitter(20e3)
const RX = GroundSampler(0:5e3:3000e3, Fields.Ez)

const DAY = Species(LMP.QE, LMP.ME, z->waitprofile(z, 75, 0.3), electroncollisionfrequency)
const NIGHT = Species(LMP.QE, LMP.ME, z->waitprofile(z, 82, 0.6), electroncollisionfrequency)
nothing  #hide

# We'll define a function to which we can pass the day or night `Species` and return the
# electric field amplitude.

function varyground(prf)
    amps = Vector{Vector{Float64}}(undef, length(GROUND))
    for i = 1:length(GROUND)
        wvg = HomogeneousWaveguide(BFIELD, prf, GROUND[i])
        E, a, p = propagate(wvg, TX, RX)
        amps[i] = a
    end
    return amps
end
nothing  #hide

# And here's a function to plot each of the curves.

function buildplots!(p, amps)
    cmap = palette(:thermal, length(GROUND)+1)  # +1 so we don't get to yellow

    for i = 1:length(GROUND)
        plot!(p, RX.distance/1000, amps[i];
              label=@sprintf("%d, %.1g", GROUND[i].ϵᵣ, GROUND[i].σ), color=cmap[i]);
    end
end
nothing  #hide

# First, the daytime results.

amps = varyground(DAY)

p = plot()
buildplots!(p, amps)
plot!(p; size=(600,400), ylims=(0, 95), title="Day", legend=(0.85, 1.02),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="ϵᵣ, σ")
#md savefig(p, "ground_day.png"); nothing # hide
#md # ![](ground_day.png)

# And now nighttime.

amps = varyground(NIGHT)

p = plot()
buildplots!(p, amps)
plot!(p; size=(600,400), ylims=(0, 95), title="Night", legend=(0.85, 1.02),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="ϵᵣ, σ")
#md savefig(p, "ground_night.png"); nothing # hide
#md # ![](ground_night.png)

# Low ground conductivity can have a significant influence on the signal propagation -
# there is strong attenuation.
# These low conductivities can be found in areas of sea or polar ice and industrial or
# city areas.
# 
# The influence of ground conductivity on the received signal has similar impact on
# the day and night scenarios.
