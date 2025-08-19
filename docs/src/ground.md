# Ground
 
LongwaveModePropagator treats the ground boundary of the Earth-ionosphere waveguide
as a homogeneous, nonmagnetic material described by its electrical conductivity ``\sigma``
and relative permittivity ``\epsilon_r``.
 
In this example we'll compare the effect of ground conductivity on amplitude along
the ground in the waveguide for day and nighttime ionospheres.

```@example ground
using Plots, Printf
using LongwaveModePropagator
const LMP = LongwaveModePropagator
nothing  # hide
```

LWPC uses standard ground indices to describe combinations of relative permittivity
and conductivity.
The standard indices can be accessed from LMP as [`GROUND`](@ref).

```@repl ground
sort(GROUND)
```

We'll define all of the constant types and the day and night profiles.

## Constant values

```@example ground
const BFIELD = BField(50e-6, π/2, 0)

const TX = Transmitter(20e3)
const RX = GroundSampler(0:5e3:3000e3, Fields.Ez)

const DAY = Species(LMP.QE, LMP.ME, z->waitprofile(z, 75, 0.3), electroncollisionfrequency)
const NIGHT = Species(LMP.QE, LMP.ME, z->waitprofile(z, 82, 0.6), electroncollisionfrequency)
nothing  # hide
```

We'll define a function to which we can pass the day or night `Species` and return the
electric field amplitude.

```@example ground
function varyground(prf)
    amps = Vector{Vector{Float64}}(undef, length(GROUND))
    for i = 1:length(GROUND)
        wvg = HomogeneousWaveguide(BFIELD, prf, GROUND[i])
        E, a, p = propagate(wvg, TX, RX)
        amps[i] = a
    end
    return amps
end
nothing  # hide
```

And here's a function to plot each of the curves.

```@example ground
function buildplots!(p, amps)
    cmap = palette(:thermal, length(GROUND)+1)  # +1 so we don't get to yellow

    for i = 1:length(GROUND)
        plot!(p, RX.distance/1000, amps[i];
              label=@sprintf("%d, %.1g", GROUND[i].ϵᵣ, GROUND[i].σ), color=cmap[i])
    end
end
nothing  # hide
```

First, the daytime results:

```@example ground
amps = varyground(DAY)

p = plot()
buildplots!(p, amps)
plot!(p; size=(600,400), ylims=(0, 95), title="Day", legend=(0.85, 1.02),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="ϵᵣ, σ")
savefig("plot-ground_a_day.svg"); nothing  # hide
```

![](plot-ground_a_day.svg)

And now nighttime:

```@example ground
amps = varyground(NIGHT)

p = plot()
buildplots!(p, amps)
plot!(p; size=(600,400), ylims=(0, 95), title="Night", legend=(0.85, 1.02),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="ϵᵣ, σ")
savefig("plot-ground_a_night.svg"); nothing  # hide
```

![](plot-ground_a_night.svg)

Low ground conductivity can have a significant influence on the signal propagation -
there is strong attenuation.
These low conductivities can be found in areas of sea or polar ice and industrial or
city areas.
 
The influence of ground conductivity on the received signal has similar impact on
the day and night scenarios.
