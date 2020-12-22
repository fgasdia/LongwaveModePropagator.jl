# # Mesh grid for mode finding - Part 2
#
# In [Mesh grid for mode finding](@ref) we created a function `trianglemesh` to
# initialize the grid used by the GRPF algorithm for mode finding. In part 2 we'll
# look at the final `grpf` grid for several different scenarios.
#
# Here are the packages needed in this example.
# Rather than redefining `trianglemesh` from part 1, here we can load a nearly
# equivalent function from `LongwaveModePropagator.jl`.

using Plots, DisplayAs
using RootsAndPoles
using VoronoiDelaunay

using ..LongwaveModePropagator
using ..LongwaveModePropagator: QE, ME, solvemodalequation, trianglemesh

## Using the GR backend
gr(legend=false);

# Here are the scenarios defined in part 1.

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

# We'll use the following triangle mesh domain.

zbl = deg2rad(complex(30.0, -10.0))
ztr = deg2rad(complex(89.9, 0.0))
Δr = deg2rad(0.5)

v = trianglemesh(zbl, ztr, Δr);

# ## Global complex roots and pole finding
#
# Let's evaluate the performance of the GRPF algorithm using `trianglemesh` to
# generate the initial mesh grid. `grpf` will adaptively refine the mesh to obtain
# a more accurate estimate of the position of the roots and poles.
#
# We will pass the `PlotData()` argument to obtain additional information on the
# function phase for plotting.

params = LMPParams().grpfparams
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, day_mid_me),
                                                  v, PlotData(), params);

# We need to do a little preprocessing before we can plot the triangle edges.

function prepareplotdata(tess, quadrants, phasediffs, g2f)
    edges = [e for e in delaunayedges(tess)]

    edgecolors = fill(5, length(edges))
    for ii in eachindex(edges)
        if phasediffs[ii] == 2  # candidate edges (but reallllly hard to see on plot)
            edgecolors[ii] = 6
        elseif phasediffs[ii] == 0
            edgecolors[ii] = quadrants[getindex(geta(edges[ii]))]
        end
    end

    edgecolors = repeat(edgecolors, inner=3)  # x, y is xₐ, xᵦ, NaN, repeat

    x, y = getplotxy(edges)
    z = g2f.(x, y)

    I = sortperm(edgecolors)

    return z[I], edgecolors[I]
end

z, edgecolors = prepareplotdata(tess, quads, phasediffs, g2f)

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
zdeg = rad2deg.(z)

pal = [
    colorant"#9E3D36",
    colorant"#C99478",
    colorant"#7599C2",
    colorant"#5C389E",
    colorant"#404040",
    RGB(0.0, 0.0, 0.0)
]

img = plot(real(zdeg), imag(zdeg), group=edgecolors, palette=pal, linewidth=1.5,
           xlims=(30, 90), ylims=(-10, 0),
           xlabel="real(θ)", ylabel="imag(θ)",
           size=(600,400),
           title=day_mid_title)
plot!(img, real(rootsdeg), imag(rootsdeg), color="red",
      seriestype=:scatter, markersize=5)
plot!(img, real(polesdeg), imag(polesdeg), color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)
DisplayAs.PNG(img)  #hide

# In the plot, roots are marked with red circles and poles are marked with red triangles.
# The automatic refinement of the mesh is clearly visible.
#
# Here are similar plots for the other three scenarios.

## Daytime ionosphere, low frequency
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, day_low_me),
                                                  v, PlotData(), params);
z, edgecolors = prepareplotdata(tess, quads, phasediffs, g2f)

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
zdeg = rad2deg.(z)

img = plot(real(zdeg), imag(zdeg), group=edgecolors, palette=pal, linewidth=1.5,
           xlims=(30, 90), ylims=(-10, 0),
           xlabel="real(θ)", ylabel="imag(θ)",
           size=(600,400),
           title=day_low_title)
plot!(img, real(rootsdeg), imag(rootsdeg), color="red",
      seriestype=:scatter, markersize=5)
plot!(img, real(polesdeg), imag(polesdeg), color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)
DisplayAs.PNG(img)  #hide

# At 100 kHz, `grpf` requires significantly more mesh refinements and takes considerably
# more time to run.

## Daytime ionosphere, high frequency
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, day_high_me),
                                                  v, PlotData(), params);
z, edgecolors = prepareplotdata(tess, quads, phasediffs, g2f)

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
zdeg = rad2deg.(z)

img = plot(real(zdeg), imag(zdeg), group=edgecolors, palette=pal, linewidth=1.5,
           xlims=(30, 90), ylims=(-10, 0),
           xlabel="real(θ)", ylabel="imag(θ)",
           size=(600,400),
           title=day_high_title)
plot!(img, real(rootsdeg), imag(rootsdeg), color="red",
      seriestype=:scatter, markersize=5)
plot!(img, real(polesdeg), imag(polesdeg), color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)
DisplayAs.PNG(img)  #hide

#

## Nighttime ionosphere, 20 kHz
roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, night_mid_me),
                                                  v, PlotData(), params);
z, edgecolors = prepareplotdata(tess, quads, phasediffs, g2f)

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
zdeg = rad2deg.(z)

img = plot(real(zdeg), imag(zdeg), group=edgecolors, palette=pal, linewidth=1.5,
           xlims=(30, 90), ylims=(-10, 0),
           xlabel="real(θ)", ylabel="imag(θ)",
           size=(600,400),
           title=night_mid_title)
plot!(img, real(rootsdeg), imag(rootsdeg), color="red",
      seriestype=:scatter, markersize=5)
plot!(img, real(polesdeg), imag(polesdeg), color="red",
      seriestype=:scatter, markershape=:utriangle, markersize=5)
DisplayAs.PNG(img)  #hide

# In `LongwaveModePropagator.jl` the [`findmodes`](@ref) function is used to build the
# initial meshgrid using `trianglemesh`, call `grpf`, and additionally run some
# uniqueness and validity checks before returning a vector of eigenangles.
