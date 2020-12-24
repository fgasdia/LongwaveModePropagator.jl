# # Wavefield integration
#
# The process of mode conversion between two different segments of `HomogeneousWaveguide`'s
# requires the computation of electromagnetic wavefields from ground to a great
# height in the ionosphere.
#
# [Pitteway](https://doi.org/10.1098/rsta.1965.0004) presents one well-known method of
# integrating the wavefields that breaks the solution for the differential wave
# equations in the anisotropic ionosphere into solutions for a _penetrating_ and
# _non-penetrating_ mode. Because of the amplitude of the two waves varies greatly
# during the integration, Pitteway scales and orthogonalizes the solutions to maintain
# their independence.
#
# In this example we'll reproduce the figures in Pitteway's paper.
#
# First we'll import the necessary packages.

using CSV, Interpolations
using Plots, DisplayAs

using ..LongwaveModePropagator
using ..LongwaveModePropagator: QE, ME, integratewavefields

# ## The ionosphere

# Pitteway uses an ionosphere presented in
# [Piggot et. al., 1965](https://doi.org/10.1098/rsta.1965.0005).
#
# ![piggott_ionosphere](Piggott_ionosphere.png)
#
# He uses the midday profile with a 16 kHz radio wave at an angle of incidence of
# 40° from normal. We'll also assume the magnetic field has a strength of 50,000 nT
# at a dip angle of 68° and azimuth of 111°.

ea = EigenAngle(deg2rad(40))
frequency = Frequency(16e3)
bfield = BField(50e-6, deg2rad(68), deg2rad(111))
ground = Ground(15, 0.001)

# To define the electron density and collision frequency profile, we will read in a
# digitized table of the curves from Piggott and linearly interpolate them.

data = CSV.File(joinpath("..", "examples", "piggott1965_data.csv"))

## interpolation object
day_itp = interpolate((data.day_ne_z,), data.day_ne, Gridded(Linear()))
day_etp = extrapolate(day_itp, Line())
dayfcn(z) = (v = day_etp(z); v > 0 ? v : 0.001)

## nu has some rows of `missing` data
filtnu_z = collect(skipmissing(data.nu_z))
filtnu = collect(skipmissing(data.nu))
nu_itp = interpolate((filtnu_z,), filtnu, Gridded(Linear()))
nu_etp = extrapolate(nu_itp, Line())
nufcn(z) = (v = nu_etp(z); v > 0 ? v : 0.001)

day = Species(QE, ME, dayfcn, nufcn)


## params = LMPParams(earthcurvature=false)

zs = 90e3:-500:0
zskm = zs/1000

e = integratewavefields(zs, ea, frequency, bfield, day, unscale=true)
e_unscaled = integratewavefields(zs, ea, frequency, bfield, day, unscale=false)

ex1 = getindex.(e, 1)
ex1_unscaled = getindex.(e_unscaled, 1)
hx2 = getindex.(e, 7)
hx2_unscaled = getindex.(e_unscaled, 7)

# The fields of e are ex1, -ey1, hx1, hy1

plot(real(ex1), zskm)
plot!(real(ex1_unscaled), zskm)

plot(real(hx2), zskm)
plot!(real(hx2_unscaled), zskm)

zs = 80e3:-500:50e3
zskm = zs/1000

e = integratewavefields(zs, ea, frequency, bfield, day, unscale=true)

ex1 = getindex.(e, 1)
ey1 = -getindex.(e, 2)
ex2 = getindex.(e, 5)
hx2 = getindex.(e, 7)

plot(real(ex1), zskm, color="black", linewidth=1.5)
plot!(-imag(ex1), zskm, color="black", linewidth=1.5, linestyle=:dash)
plot!(abs.(ex1), zskm, color="black", linewidth=3)
plot!(-abs.(ex1), zskm, color="black", linewidth=3)

plot(real(ey1), zskm, color="black", linewidth=1.5)
plot!(-imag(ey1), zskm, color="black", linewidth=1.5, linestyle=:dash)
plot!(abs.(ey1), zskm, color="black", linewidth=3)
plot!(-abs.(ey1), zskm, color="black", linewidth=3)
