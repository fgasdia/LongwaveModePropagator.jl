using Random, Dates, Statistics
using Plots
using Plots.Measures
using RootsAndPoles, Interpolations
using ProgressMeter

using GeographicLib
using LMPTools
import FaradayInternationalReferenceIonosphere as FIRI

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

###
# Generate random segmented scenarios

function randscenarios(N, freq, dt, Δr)
      rng = MersenneTwister(1234)
      x = 0:200e3:2000e3
      dists = 0:1e3:2000e3

      # Mesh
      zbl = complex(deg2rad(10), deg2rad(-30))
      ztr = complex(deg2rad(89.9), 0.0)
      mesh = trianglemesh(zbl, ztr, Δr)

      Nwvg = SegmentedWaveguide[]
      Namp = Matrix{Float64}(undef, 2001, N)
      Nphase = similar(Namp)
      for n = 1:N
            # Random tx position between 0 and 60°N.
            txlat = 59*rand(rng).+1
            txlon = 359*rand(rng)
            tx = Transmitter("tx", txlat, txlon, freq)

            n == 1 && @assert txlat == 35.8598336812769  # because of rng

            # Random bearing
            geoaz = 360*rand(rng)
            line = GeodesicLine(txlat, txlon; azi=geoaz)

            # Actual IGRF, ground map, FIRI profile
            wvgs = Vector{HomogeneousWaveguide{Species}}(undef, 11)
            for i = 1:11
                  wpt = forward(line, x[i])
                  bfield = igrf(geoaz, wpt.lat, wpt.lon, year(dt))
                  sza = zenithangle(wpt.lat, wpt.lon, dt)
                  ground = GROUND[get_groundcode(wpt.lat, wpt.lon)]
                  elprofile = FIRI.extrapolate(FIRI.firi(sza, abs(wpt.lat)), 0:1e3:110e3)
                  ne = Interpolations.interpolate(0:1e3:110e3, elprofile, FritschButlandMonotonicInterpolation())
                  species = Species(QE, ME, ne, electroncollisionfrequency)
                  wvgs[i] = HomogeneousWaveguide(bfield, species, ground, x[i])
            end
            wvg = SegmentedWaveguide(wvgs)
            push!(Nwvg, wvg)

            E, amplitude_db, phase_rad = propagate(wvg, tx, GroundSampler(dists, Fields.Ez); mesh)
            Namp[:,n] .= amplitude_db
            Nphase[:,n] .= rad2deg.(phase_rad)
      end

      return Nwvg, Namp, Nphase
end

dt = DateTime(2022, 7, 7)
freqs = vcat(3e3:1e3:30e3, 32e3:2e3:60e3, 65e3:5e3:100e3)
Δrs = deg2rad.([2, 1, 0.5, 0.4, 0.3, 0.2, 0.1])
# freqs = (5e3, 6e3, 7e3)
# Δrs = deg2rad.([3, 2, 1])
wvgs = Matrix{Vector{SegmentedWaveguide}}(undef, length(Δrs), length(freqs))
amps = Matrix{Matrix{Float64}}(undef, length(Δrs), length(freqs))
phases = similar(amps)
@showprogress for j in eachindex(freqs)
      for i in eachindex(Δrs)
            Nwvg, Namp, Nphase = randscenarios(5, freqs[j], dt, Δrs[i]);
            wvgs[i,j] = Nwvg
            amps[i,j] = Namp
            phases[i,j] = Nphase
      end
end
damps = amps .- reshape(amps[end,:], 1, length(freqs))
netda = mean.(abs, damps; dims=1)

cmap = cgrad(:matter, length(Δrs)-1, categorical = true)
scatter()
for i = 1:length(Δrs)-1
      for j in eachindex(freqs)
            if j == 1
                  scatter!(fill(freqs[j], 5), netda[i,j]', color=cmap[i], label=round(rad2deg(Δrs[i]), digits=1))
            else
                  scatter!(fill(freqs[j], 5), netda[i,j]', color=cmap[i], label=nothing)
            end
      end
end
scatter!(label=nothing)
