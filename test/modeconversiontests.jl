using Test
using LinearAlgebra
using Statistics
using StaticArrays
# using Plots
# using NumericalIntegration
# using Trapz  # for testing only
using Parameters
using CSV
using DataFrames

using RootsAndPoles

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C


# function test_calculate_wavefields!()
#     bfield, tx, ground, electrons, ea, zs = resonant_scenario()
#
#     wavefields = Wavefields(ea, zs)
#     adjoint_wavefields = Wavefields(ea, zs)
#
#     calculate_wavefields!(wavefields, adjoint_wavefields,
#                           bfield, tx.frequency, ground, electrons)
#
#     return wavefields, adjoint_wavefields
# end


function homoscenario(threshold)
    waveguide =  HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 75, 0.32, threshold=threshold),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001))

    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    rx = GroundSampler(0:5e3:2000e3, LWMS.FC_Ez)

    E, phase, amp = LWMS.bpm(waveguide, tx, rx)

    return E, phase, amp
end


function mc_scenario(threshold)
    waveguide = LWMS.SegmentedWaveguide(HomogeneousWaveguide)

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 75, 0.32, threshold=threshold),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001)))

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 70, 0.25, threshold=threshold),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001), 1000e3))


    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    rx = GroundSampler(0:5e3:2000e3, LWMS.FC_Ez)

    E, phase, amp = LWMS.bpm(waveguide, tx, rx)

    return E, phase, amp
end




# i = 6
# EH = reshape(reinterpret(ComplexF64, EH), 6, :)
# plot(real(EH[i,:]), zs/1000, color="black")
# plot!(imag(EH[i,:]), zs/1000, color="black")
# plot!(abs.(EH[i,:]), zs/1000, color="black")
# plot!(-abs.(EH[i,:]), zs/1000, color="black")
#
# EHadjoint = reshape(reinterpret(ComplexF64, EHadjoint), 6, :)
# plot!(real(EHadjoint[i,:]), zs/1000, color="blue")
# plot!(imag(EHadjoint[i,:]), zs/1000, color="blue")
# plot!(abs.(EHadjoint[i,:]), zs/1000, color="blue")
# plot!(-abs.(EHadjoint[i,:]), zs/1000, color="blue")
#
#
# x = 0.28656210625768974
# @test LWMS.pow23(x) == LWMS.pow23(complex(x))
# y = 100*x
# @test LWMS.pow23(y) == LWMS.pow23(complex(y))




using GRUtils
colorscheme("dark")  # matches atom, otherwise "light"

ceilto(x, a) = ceil(Int, x/a)*a
floorto(x, a) = floor(Int, x/a)*a

basepath = "/home/forrest/research/LAIR/ModeSolver/lwpc_comparisons/"

tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
rx = GroundSampler(0:5e3:2000e3, LWMS.FC_Ez)
X = LWMS.distance(rx,tx);

# Homogeneous scenario
E, phase, amp = homoscenario(3e9)
namp = copy(amp)

# TEMP  first point in amp is NaN which causes it to not plot
namp[1] = namp[2]

rawdft = CSV.File(joinpath(basepath, "homogeneousionosphere", "dfts.csv");
                  header=false)
dft = DataFrame(dist=rawdft.Column1, amp=rawdft.Column2, phase=rawdft.Column3);
dft = dft[1:findfirst(dft.dist .== X[end]/1000),:]
dft = dft[1:10:end,:]
dft.amp .+= 114  # used for Jamesina comparison


subplot(4, 1, (1, 3));
plot(X/1000, namp, linewidth=1.5)
oplot(X/1000, dft.amp, linewidth=1.5)
legend("BPM", "FDTD");
ylabel("Amp (dB μV/m)");

dftdiff = namp-dft.amp
ylim_min = floor(quantile(dftdiff, 0.01))
ylim_max = ceil(quantile(dftdiff, 0.99))

lowplot = subplot(4, 1, 4);
plot(X/1000, dftdiff, "r-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
ylabel("Δ BPM");
xlabel("Distance (km)");
gcf()



# Single transition
scenario = "singletransition"
raw = CSV.File(joinpath(basepath, scenario, "singletransition.log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]

rawdft = CSV.File(joinpath(basepath, scenario, "dfts.csv");
                  header=false)
dft = DataFrame(dist=rawdft.Column1, amp=rawdft.Column2, phase=rawdft.Column3);
dft = dft[1:findfirst(dft.dist .== X[end]/1000),:]
dft = dft[1:10:end,:]
dft.amp .+= 114  # used for Jamesina comparison


_, thresh_phase, thresh_amp = mc_scenario(3e9)
nthresh_amp = copy(thresh_amp)
nthresh_amp[1] = nthresh_amp[2]

subplot(4, 1, (1, 3));
plot(X/1000, nthresh_amp, linewidth=1.5);
oplot(X/1000, dat.amp, linewidth=1.5);
oplot(X/1000, dft.amp, linewidth=1.5);
legend("BPM_thresh", "LWPC", "FDTD");
ylabel("Amp (dB μV/m)");

dftdiff = nthresh_amp .- dft.amp
datdiff = nthresh_amp .- dat.amp
ylim_min = floorto(min(quantile(dftdiff, 0.01), quantile(datdiff, 0.01)), 5) - 1
ylim_max = ceilto(max(quantile(dftdiff, 0.99), quantile(datdiff, 0.99)), 5) + 1

lowplot = subplot(4, 1, 4);
plot(X/1000, datdiff, "r-", linewidth=1.5);
oplot(X/1000, dftdiff, "y-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
ylabel("Δ BPM");
xlabel("Distance (km)");
gcf()


ylim_min = min(quantile(rad2deg.(thresh_phase), 0.01), quantile(dat.phase, 0.01), quantile(dft.phase, 0.01))
ylim_max = max(quantile(rad2deg.(thresh_phase), 0.99), quantile(dat.phase, 0.99), quantile(dft.phase, 0.99))

ylim_min = floorto(ylim_min, 5) - 1
ylim_max = ceilto(ylim_max, 5) + 1

subplot(4, 1, (1, 3));
plot(X/1000, rad2deg.(thresh_phase), linewidth=1.5);
oplot(X/1000, dat.phase, linewidth=1.5);
oplot(X/1000, dft.phase, linewidth=1.5);
ylim(ylim_min, ylim_max)
legend("BPM_thresh", "LWPC", "FDTD");
ylabel("Phase (deg)");

dftdiff = rad2deg.(thresh_phase)-dft.phase
datdiff = rad2deg.(thresh_phase)-dat.phase
ylim_min = floorto(min(quantile(dftdiff, 0.01), quantile(datdiff, 0.01)), 5) - 1
ylim_max = ceilto(max(quantile(dftdiff, 0.99), quantile(datdiff, 0.99)), 5) + 1

subplot(4, 1, 4);
plot(dat.dist, datdiff, "r-", linewidth=1.5);
oplot(dat.dist, dftdiff, "y-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
xlabel("Distance (km)");
ylabel("Δ BPM");
gcf()
