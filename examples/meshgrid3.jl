using Plots
using Plots.Measures
using RootsAndPoles
using ProgressMeter

using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, solvemodalequation, trianglemesh

function modeequationphase(me, mesh)
      phase = Vector{Float64}(undef, length(mesh))
      Threads.@threads for i in eachindex(mesh)
            f = solvemodalequation(deg2rad(mesh[i]), me)
            phase[i] = rad2deg(angle(f))
      end
      return phase
end

twilightquads = [
    colorant"#9E3D36",
    colorant"#C99478",
    colorant"#7599C2",
    colorant"#5C389E",
    colorant"#404040",
    RGB(0.0, 0.0, 0.0)
]

elowfrequency = Frequency(3e3)
lowfrequency = Frequency(8e3)
midfrequency = Frequency(20e3)
highfrequency = Frequency(100e3)

day = Species(QE, ME, z->waitprofile(z, 70, 0.25), electroncollisionfrequency)
night = Species(QE, ME, z->waitprofile(z, 85, 0.8), electroncollisionfrequency)

daywaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), day, GROUND[5])
nightwaveguide = HomogeneousWaveguide(BField(50e-6, π/2, 0), night, GROUND[5])

day_elow_me = PhysicalModeEquation(elowfrequency, daywaveguide)
day_low_me = PhysicalModeEquation(lowfrequency, daywaveguide)
day_mid_me = PhysicalModeEquation(midfrequency, daywaveguide)
day_high_me = PhysicalModeEquation(highfrequency, daywaveguide)

night_low_me = PhysicalModeEquation(lowfrequency, nightwaveguide)
night_mid_me = PhysicalModeEquation(midfrequency, nightwaveguide)
night_high_me = PhysicalModeEquation(highfrequency, nightwaveguide)

day_elow_title = "3 kHz\nh′: 75, β: 0.35"
day_low_title = "80 kHz\nh′: 75, β: 0.35"
day_mid_title = "20 kHz\nh′: 75, β: 0.35"
day_high_title = "100 kHz\nh′: 75, β: 0.35"
night_mid_title = "20 kHz\nh′: 85, β: 0.9"


Δr = 0.2
x = 10:Δr:89.9
y = -30:Δr:0
grid = x .+ 1im*y';

scn = night_mid_me
scn_title = night_mid_title

params = LMPParams().grpfparams

mesh = rectangulardomain(complex(first(x), first(y)), complex(last(x), last(y)), Δr)


phase = modeequationphase(scn, grid);

img = heatmap(x, y, reshape(phase, length(x), length(y))';
      color=:twilight, clims=(-180, 180),
      xlims=(10, 90), ylims=(-30, 0),
      xlabel="real(θ)", ylabel="imag(θ)",
      title=scn_title, label=nothing,
      right_margin=2mm)


roots, poles, quads, phasediffs, tess, g2f = grpf(θ->solvemodalequation(θ, scn),
                                                  deg2rad.(mesh), PlotData(), params);

z, edgecolors = getplotdata(tess, quads, phasediffs, g2f)

rootsdeg = rad2deg.(roots)
polesdeg = rad2deg.(poles)
zdeg = rad2deg.(z)

# img = plot(real(zdeg), imag(zdeg); group=edgecolors, palette=twilightquads, linewidth=1.5,
#            xlims=(30, 90), ylims=(-10, 0),
#            xlabel="real(θ)", ylabel="imag(θ)", legend=false,
#            title=day_mid_title);
plot!(img, real(rootsdeg), imag(rootsdeg); color="red",
      seriestype=:scatter, markersize=5, label=nothing)
# plot!(img, real(polesdeg), imag(polesdeg); color="red",
#       seriestype=:scatter, markershape=:utriangle, markersize=5)


function rootextent(freqs, wvg, mesh)
      ext = Vector{ComplexF64}(undef, length(freqs))
      @showprogress for i in eachindex(freqs)
            me = PhysicalModeEquation(Frequency(freqs[i]), wvg)
            roots, _ = grpf(θ->solvemodalequation(θ, me), deg2rad.(mesh), params)
            rootsdeg = rad2deg.(roots)
            ext[i] = complex(minimum(real, rootsdeg), minimum(imag, rootsdeg))
      end
      return ext
end
freqs = vcat(3e3:1e3:30e3, 32e3:2e3:60e3, 65e3:5e3:100e3)
exts = rootextent(freqs, daywaveguide, mesh)
nexts = rootextent(freqs, nightwaveguide, mesh)
scatter(real(exts), imag(exts), zcolor=freqs/1e3)
