using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
const LMP = LongwaveModePropagator


f = 50e3
tx = Transmitter(f)

ranges = 0:5e3:2000e3
field = Fields.Ez
rx = GroundSampler(ranges, field)

bfield = BField(50e-6, π/2, 0)
electrons = Species(QE, ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)
ground = GROUND[5]

waveguide = HomogeneousWaveguide(bfield, electrons, ground)
me = PhysicalModeEquation(tx.frequency, waveguide)


mesh = deg2rad.([0-40.0im, 0-0.0im, 89.99-40.0im, 89.99-0.0im])
modes = findmodes(me, mesh)
E1, a1, p1 = propagate(waveguide, tx, rx; modes=modes)


mesh = deg2rad.([0-30.0im, 0-0.0im, 89.99-30.0im, 89.99-0.0im])
modes = findmodes(me, mesh)
E2, a2, p2 = propagate(waveguide, tx, rx; modes=modes)


mesh = deg2rad.([0-20.0im, 0-0.0im, 89.99-20.0im, 89.99-0.0im])
modes = findmodes(me, mesh)
E3, a3, p3 = propagate(waveguide, tx, rx; modes=modes)


all(<(0.1), a3[2:end] - a1[2:end])
all(<(0.1), p3[2:end] - p1[2:end])

using Plots

Plots.plot(ranges/1000, [a1 a2 a3];
     xlabel="range (km)", ylabel="amplitude (dB)",
     linewidth=1.5)

Plots.plot(ranges/1000, [p1 p2 p3];
     xlabel="range (km)", ylabel="phase (rad)",
     linewidth=1.5)



using Contour

reals = deg2rad.(0:89)
imags = deg2rad.(-60:0)

freqs = [5e3, 10e3, 20e3, 40e3, 80e3, 120e3]
maxi = zeros(Float64, length(freqs))
for (i, f) in enumerate(freqs)
    z = [attenuation(complex(x, y), Frequency(f)) for x in reals, y in imags]
    cl = Contour.contour(reals, imags, z, 1000)
    coords = Contour.coordinates(only(Contour.lines(cl)))
    maxi[i] = maximum(coords[2])
end

# solving attenuation for max imag theta with atten = 10^3
neper2dB = 20log10(exp(1))
rad2deg.(asinh.(10^(-3) ./ (neper2dB.*ks)))


Plots.contourf(reals, imags, (x, y)->log10(attenuation(complex(x, y), tx.frequency));
     xlabel="real(θ) (rad)", ylabel="imag(θ) (rad)",
     colorbar_title="log10(attenuation) (dB)",
     linewidth=1)

Plots.plot!(coords, color="white", linestyle=:dash)
