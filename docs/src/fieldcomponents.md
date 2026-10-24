# Multiple field components

Multiple components of the propagating electromagnetic field can be calculated and sampled
at once. This is greatly advantageous compared to calling `propagate` multiple times
because the runtime of `propagate` is dominated (>90%) by searching for waveguide
eigenangles. Relatively little compute time is required to calculate and sum the fields
associated with each individual waveguide mode.
 
Let's load the necessary packages.

```@example components
using Plots

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
```

## [`Fields`](@ref)

Each supported electromagnetic field is defined in an `Enum` called `Fields`. The supported fields are:

```@repl components
Fields.Field
```

The ``x`` axis extends from the transmitter along the ground in the direction of the
receiver. The ``z`` axis extends vertically upward into the ionosphere so that the
wavefields are propagating in the ``x-z`` plane. The ``y`` axis is perpendicular to the
plane and completes the right handed coordinate system.

Individual fields are specified with a dot syntax:

```@repl components
Fields.Ez
```

`Fields.Ez`, `Fields.Ex`, and `Fields.Ey` are self explanatory. `Fields.E` returns all
three electric field components at once.

!!! warning

     It is not recommended to rely upon the integer representation of the `Enum`.

## Propagating multiple fields: HomogeneousWaveguide

Using the same conditions defined in [the introduction](@ref basic_homogeneous), we'll
define a model using a new `GroundSampler` with fieldcomponent `Fields.E`, run the model,
and plot the results.

```@example components
f = 24e3
tx = Transmitter(f)

ranges = 0:1e3:2000e3
field = Fields.E
rx = GroundSampler(ranges, field)

bfield = BField(50e-6, π/2, 0)
electrons = Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)
ground = GROUND[5]

waveguide = HomogeneousWaveguide(bfield, electrons, ground)
```

The [`propagate`](@ref) function returns a tuple of complex electric field,
amplitude in dB μV/m, and phase in radians.

```@example components
E, a, p = propagate(waveguide, tx, rx)
nothing  # hide
```

You'll notice that `E`, `a`, and `p` are 2001 × 3 matrices

```@example components
size(E), size(a), size(p)
```

The columns are, in order, the `Ez`, `Ey`, and `Ex` field components.

Here are quick plots of the amplitude:

```@example components
fieldlabels = ["Ez" "Ey" "Ex"]

plot(ranges/1000, a;
     xlabel="range (km)", ylabel="amplitude (dB)",
     linewidth=1.5, label=fieldlabels)
savefig("plot-fieldcomponents_a.svg"); nothing  # hide
```

![](plot-fieldcomponents_a.svg)

and phase:

```@example components
plot(ranges/1000, rad2deg.(p);
     xlabel="range (km)", ylabel="phase (deg)",
     linewidth=1.5, label=fieldlabels)
savefig("plot-fieldcomponents_p.svg"); nothing  # hide
```

![](plot-fieldcomponents_p.svg)

The `Ey` phase grows rapidly - an alternative plot would `mod2pi` the results to
"undo" the phase unwrapping applied by `propagate`.

```@example components
plot(ranges/1000, rad2deg.(mod2pi.(p));
     xlabel="range (km)", ylabel="phase (deg)",
     linewidth=1.5, label=fieldlabels)
savefig("plot-fieldcomponents_p_mod.svg"); nothing  # hide
```

![](plot-fieldcomponents_p_mod.svg)

## Propagating multiple fields: SegmentedWaveguide

In this example we'll have a [`SegmentedWaveguide`](@ref) with two segments.

```@example components
distances = [0.0, 1000e3]
species = [Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency),
           Species(QE, ME, z->waitprofile(z, 82, 0.5), electroncollisionfrequency)]

waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield, species[i], ground,
                                                     distances[i]) for i in 1:2])
```

We can [`propagate`](@ref) just as before

```@example components
E, a, p = propagate(waveguide, tx, rx)
nothing  # hide
```

Here are quick plots of the amplitude:

```@example components
plot(ranges/1000, a;
     xlabel="range (km)", ylabel="amplitude (dB)",
     linewidth=1.5, label=fieldlabels)
savefig("plot-fieldcomponents_a_segm.svg"); nothing  # hide
```

![](plot-fieldcomponents_a_segm.svg)

and phase:

```@example components
plot(ranges/1000, rad2deg.(p);
     xlabel="range (km)", ylabel="phase (deg)",
     linewidth=1.5, label=fieldlabels)
savefig("plot-fieldcomponents_p_segm.svg"); nothing  # hide
```

![](plot-fieldcomponents_p_segm.svg)

and again with phase wrapping:

```@example components
plot(ranges/1000, rad2deg.(mod2pi.(p));
     xlabel="range (km)", ylabel="phase (deg)",
     linewidth=1.5, label=fieldlabels)
savefig("plot-fieldcomponents_p_mod_segm.svg"); nothing  # hide
```

![](plot-fieldcomponents_p_mod_segm.svg)
