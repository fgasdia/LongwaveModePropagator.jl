# Mesh grid for mode finding

The [Global complex Roots and Poles Finding](https://doi.org/10.1109/TAP.2018.2869213)
(GRPF) algorithm searches for roots and poles of complex-valued functions by sampling
the function at the nodes of a regular mesh and analyzing the function's complex phase.
The mesh search region can be arbitrarily shaped and guarantees that none of
the roots/poles can be missed as long as the mesh step is sufficiently small.
Specifically, the phase change of the function value cannot exceed three quadrants
across mesh edges. Effectively this means there cannot be equal numbers of roots and
poles within the boundary of each identified candidate region.

Although the GRPF algorithm uses an adaptive meshing process, its ability to find roots and poles was sensitive to the initial mesh node spacing and there is no established method for _a priori_ estimation of the required initial sampling for an arbitrary function. The [Self-Adaptive Mesh Generator for Global Complex Roots and Poles Finding Algorithm](https://ieeexplore.ieee.org/document/10025853/) (SA-GRPF) can start from just four corner points on the complex plane, using a gradient calculation to identify areas that need mesh refinement. The SA-GRPF algorithm is implemented in [RootsAndPoles.jl](https://github.com/fgasdia/RootsAndPoles.jl), and is used by LongwaveModePropagator v0.5 to identify eigenangles of [`HomogeneousWaveguide`](@ref)'s.

## Mode equation

The criteria for mode propagation in the Earth-ionosphere waveguide is described by the
fundamental equation of mode theory

```math
\det(\overline{\bm{R}}(\theta)\bm{R}(\theta) - \bm{I}) = 0
```

where ``\overline{\bm{R}}(\theta)`` is the reflection coefficient matrix for the
ground and ``\bm{R}(\theta)`` is the reflection coefficient of the ionosphere.
Both are functions of the complex angle ``\theta``. The discrete values of
``\theta`` for which the fundamental equation is satisfied are known as
_eigenangles_.

First, import the necessary packages.

```@example meshgrid
using Plots
using Plots.Measures
using RootsAndPoles
const RP = RootsAndPoles

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, solvemodalequation, trianglemesh
```

## Exploring the mode equation phase

The U.S. Navy's Long Wavelength Propagation Capability (LWPC) searches the region of the
complex plane from 30° to 90° in the real axis and 0° to -10° in the
imaginary axis. The lowest order modes are nearest 90° - i0°, excluding 90°.
Modes with large negative imaginary components are highly attenuated and have
relatively little affect on the total field in the waveguide.

Let's begin by exploring the modal equation phase on a fine grid across the bottom
right complex quadrant.

We'll define five different [`PhysicalModeEquation`](@ref)'s that we'll use
throughout this example.

```@example meshgrid
verylowfrequency = Frequency(1e3)
lowfrequency = Frequency(10e3)
midfrequency = Frequency(20e3)
highfrequency = Frequency(100e3)

day = Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)
night = Species(QE, ME, z->waitprofile(z, 85, 0.9), electroncollisionfrequency)

daywaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), day, GROUND[5])
nightwaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), night, GROUND[5])

day_verylow_me = PhysicalModeEquation(verylowfrequency, daywaveguide)
day_low_me = PhysicalModeEquation(lowfrequency, daywaveguide)
day_mid_me = PhysicalModeEquation(midfrequency, daywaveguide)
day_high_me = PhysicalModeEquation(highfrequency, daywaveguide)
night_mid_me = PhysicalModeEquation(midfrequency, nightwaveguide)

day_verylow_title = "1 kHz\nh′: 75, β: 0.35"
day_low_title = "10 kHz\nh′: 75, β: 0.35"
day_mid_title = "20 kHz\nh′: 75, β: 0.35"
day_high_title = "100 kHz\nh′: 75, β: 0.35"
night_mid_title = "20 kHz\nh′: 85, β: 0.9"
nothing  # hide
```

Here we define a dense rectangular mesh.

```@example meshgrid
Δr = 0.2
x = 0:Δr:90
y = -40:Δr:0
mesh = x .+ im*y'
nothing  # hide
```

Now we simply iterate over each node of the mesh, evaluating the modal equation with
[`LongwaveModePropagator.solvemodalequation`](@ref), explicitly imported from
LongwaveModePropagator above.
We use Julia's `Threads.@threads` multithreading capability to speed up the computation.

```@example meshgrid
function modeequationphase(me, mesh)
    phase = Vector{Float64}(undef, length(mesh))
    Threads.@threads for i in eachindex(mesh)
        f = solvemodalequation(deg2rad(mesh[i]), me)
        phase[i] = rad2deg(angle(f))
    end
    return phase
end

phase = modeequationphase(day_mid_me, mesh)
nothing  # hide
```

If an instability occurs in the integration of the reflection coefficient,
it likely happened at the angle 90° + i0°, which is not a valid mode.
This isn't a problem for plotting though!

Plotting the results, we can see that there are clearly identifiable locations where
white, black, blue, and orange, each representing a different quadrant of the
complex plane, all meet.
Each of these locations are either a root or pole in the daytime ionosphere.

```@example meshgrid
heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-40, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_mid_title,
        right_margin=2mm)

savefig("plot-meshgrid_day_mid_me.png"); nothing  # hide
```

![](plot-meshgrid_day_mid_me.png)

We can zoom in to the upper right corner of the plot to see the lowest order modes:

```@example meshgrid
heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(30, 90), ylims=(-10, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_mid_title,
        right_margin=2mm)

savefig("plot-meshgrid_day_mid_me_zoom.png"); nothing  # hide
```

![](plot-meshgrid_day_mid_me_zoom.png)

If we switch to a nighttime ionosphere with a high Wait β parameter, we see that the
roots move closer to the axes.
A perfectly reflecting conductor has eigenangles along the real and complex axes.

```@example meshgrid
phase = modeequationphase(night_mid_me, mesh);

heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-40, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=night_mid_title,
        right_margin=2mm)

savefig("plot-meshgrid_night_mid_me.png"); nothing  # hide
```

![](plot-meshgrid_night_mid_me.png)

At lower frequencies, the roots/poles move further apart.

```@example meshgrid
phase = modeequationphase(day_low_me, mesh);

heatmap(x, y, reshape(phase, length(x), length(y))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-40, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_low_title,
        right_margin=2mm)

savefig("plot-meshgrid_day_low_me.png"); nothing  # hide
```

![](plot-meshgrid_day_low_me.png)

At even lower frequencies into ELF, we need to expand the domain of the mesh grid.

```@example meshgrid
y_verylow = -90:Δr:0
mesh_verylow = x .+ im*y_verylow';

phase = modeequationphase(day_verylow_me, mesh_verylow);

heatmap(x, y_verylow, reshape(phase, length(x), length(y_verylow))';
        color=:twilight, clims=(-180, 180),
        xlims=(0, 90), ylims=(-90, 0),
        xlabel="real(θ)", ylabel="imag(θ)",
        title=day_verylow_title,
        right_margin=2mm)

savefig("plot-meshgrid_day_verylow_me.png"); nothing  # hide
```

![](plot-meshgrid_day_verylow_me.png)

## Global complex roots and poles finding

Let's apply `rootsandpoles` to the modal equation.

[`LMPParams`](@ref) is a struct use for passing parameters across LongwaveModePropagator.
Default values are automatically inserted when instantiating the struct.
We'll extract the `finderparams` field which contains a `FinderParams` struct for passing
to `rootsandpoles`.

```@example meshgrid
initmesh = deg2rad.([0-40.0im, 0-0.0im, 89.99-40.0im, 89.99-0.0im])
mesh = ComplexMesh(initmesh)

params = FinderParams(tol=1e-5, numtasks=Threads.nthreads())
roots, poles = rootsandpoles(θ->solvemodalequation(θ, day_mid_me), mesh; params);
```

We'll write a function to help us draw each edge of the mesh.

!!! warning

    The `ComplexMesh` `mesh` is updated in place by `rootsandpoles`.

```@example meshgrid
function formatquadrants(mesh)
    xs = Vector{Float64}()
    ys = similar(xs)
    cs = Vector{Int}()
    for e in RP.each_solid_edge(mesh)
        u, v = RP.DT.edge_vertices(e)
        x, y = RP.coord(mesh, u)
        p, q = RP.coord(mesh, v)

        qu, qv = RP.quadrant(mesh, u, v)
        ΔQ = abs(qu - qv)
        ΔQ == 2 ? c = 5 :
            ΔQ == 0 ? c = RP.quadrant(mesh, u) :
                c = 0

        push!(xs, x, p, NaN)
        push!(ys, y, q, NaN)
        push!(cs, c, c, c)
    end
    return xs, ys, cs
end

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)

xs, ys, cs = formatquadrants(mesh)

twilightquads = [
    colorant"#404040",  # bndry
    colorant"#9E3D36",
    colorant"#C99478",
    colorant"#7599C2",
    colorant"#5C389E",
    RGB(0.0, 0.0, 0.0)  # candidate
]

img = plot(rad2deg.(xs), rad2deg.(ys); group=cs, palette=twilightquads, linewidth=1.5,
           xlims=(0, 90), ylims=(-40, 0),
           xlabel="real(θ)", ylabel="imag(θ)", legend=false,
           title=day_mid_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
      seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)

savefig("plot-meshgrid_day_mid_mesh.svg"); nothing  # hide
```

![](plot-meshgrid_day_mid_mesh.svg)

In the plot above, roots are marked with red circles and poles are marked
with red triangles. The mesh was automatically refined near roots and poles.

Here are similar plots for the other scenarios.

## Daytime ionosphere, 10 kHz

```@example meshgrid
mesh = ComplexMesh(initmesh)
roots, poles = rootsandpoles(θ->solvemodalequation(θ, day_low_me), mesh; params);

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
xs, ys, cs = formatquadrants(mesh)

img = plot(rad2deg.(xs), rad2deg.(ys); group=cs, palette=twilightquads, linewidth=1.5,
           xlims=(0, 90), ylims=(-40, 0),
           xlabel="real(θ)", ylabel="imag(θ)", legend=false,
           title=day_low_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
      seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)

savefig("plot-meshgrid_day_low_mesh.svg"); nothing  # hide
```

![](plot-meshgrid_day_low_mesh.svg)

## Nighttime ionosphere, 20 kHz

```@example meshgrid
mesh = ComplexMesh(initmesh)
roots, poles = rootsandpoles(θ->solvemodalequation(θ, night_mid_me), mesh; params);

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
xs, ys, cs = formatquadrants(mesh)

img = plot(rad2deg.(xs), rad2deg.(ys); group=cs, palette=twilightquads, linewidth=1.5,
           xlims=(0, 90), ylims=(-40, 0),
           xlabel="real(θ)", ylabel="imag(θ)", legend=false,
           title=night_mid_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
      seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)

savefig("plot-meshgrid_night_mid_mesh.svg"); nothing  # hide
```

![](plot-meshgrid_night_mid_mesh.svg)

At 100 kHz, `grpf` requires more mesh refinements and takes longer to run.

## Daytime ionosphere, 100 kHz

```@example meshgrid
mesh = ComplexMesh(initmesh)
roots, poles = rootsandpoles(θ->solvemodalequation(θ, day_high_me), mesh; params);

@info "$(length(roots)) roots and $(length(poles)) poles"

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
xs, ys, cs = formatquadrants(mesh)

img = plot(rad2deg.(xs), rad2deg.(ys); group=cs, palette=twilightquads, linewidth=1.5,
           xlims=(0, 90), ylims=(-40, 0),
           xlabel="real(θ)", ylabel="imag(θ)", legend=false,
           title=day_high_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
      seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)

savefig("plot-meshgrid_day_high_mesh.svg"); nothing  # hide
```

![](plot-meshgrid_day_high_mesh.svg)

Increasing the number of adaptive and total mesh nodes results in a finding several extra roots.

```@example meshgrid
newparams = FinderParams(tol=1e-5, maxadaptivenodes=10000, maxnodes=40000, numtasks=Threads.nthreads())
roots, poles = rootsandpoles(θ->solvemodalequation(θ, day_high_me), mesh; params=newparams);

@info "$(length(roots)) roots and $(length(poles)) poles"
```

## Daytime ionosphere, 1 kHz

The mesh grid does not need to be as dense at lower VLF and ELF, where the roots and poles are widely separated.

Below 8 kHz, the full lower right quadrant of the complex plane is searched for modes.
Consider additionally providing your own `params` struct with a lower value for `maxadaptivenodes` if regularly modeling lower VLF or ELF propagation.

```@example meshgrid
initmesh = deg2rad.([0-90.0im, 0-0.0im, 89.99-90.0im, 89.99-0.0im]) 
mesh = ComplexMesh(initmesh)
roots, poles = rootsandpoles(θ->solvemodalequation(θ, day_verylow_me), mesh; params);

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
xs, ys, cs = formatquadrants(mesh)

img = plot(rad2deg.(xs), rad2deg.(ys); group=cs, palette=twilightquads, linewidth=1.5,
           xlims=(0, 90), ylims=(-90, 0),
           xlabel="real(θ)", ylabel="imag(θ)", legend=false,
           title=day_verylow_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
      seriestype=:scatter, markersize=5);
plot!(img, real(polesdeg), imag(polesdeg); color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)

savefig("plot-meshgrid_day_verylow_mesh.svg"); nothing  # hide
```

![](plot-meshgrid_day_verylow_mesh.svg)
