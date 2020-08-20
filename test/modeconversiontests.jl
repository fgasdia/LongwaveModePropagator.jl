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


function homoscenario2(threshold)
    waveguide =  HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 70, 0.25, threshold=threshold),
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

function mc_scenario_3(threshold)
    waveguide = LWMS.SegmentedWaveguide(HomogeneousWaveguide)

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 75, 0.32, threshold=threshold),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001)))

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 73, 0.30, threshold=threshold),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001), 800e3))

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



using GRUtils
colorscheme("light")  # matches atom, otherwise "light"

ceilto(x, a) = ceil(Int, x/a)*a
floorto(x, a) = floor(Int, x/a)*a

basepath = "/home/forrest/research/LAIR/ModeSolver/lwpc_comparisons/"

tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
rx = GroundSampler(0:5e3:2000e3, LWMS.FC_Ez)
X = LWMS.distance(rx,tx);


# Homogeneous scenario
scenario = "homoscenario"

E, phase, amp = homoscenario(3e9)
namp = copy(amp)

# TEMP  first point in amp is NaN which causes it to not plot
namp[1] = namp[2]

rawdft = CSV.File(joinpath(basepath, scenario, "alldfts.csv");
                  header=false)
dft = DataFrame(dist=rawdft.Column1, amp=rawdft.Column2, phase=rawdft.Column3);
dft = dft[1:findfirst(dft.dist .== X[end]/1000),:]
dft = dft[1:10:end,:]
dft.amp .+= 116  # used for Jamesina comparison
dft.phase .-= 33

raw = CSV.File(joinpath(basepath, scenario, scenario*".log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]

subplot(4, 1, (1, 3));
plot(X/1000, namp, "b", linewidth=1.5)
oplot(X/1000, dat.amp, "r", linewidth=1.5)
oplot(X/1000, dft.amp, "y", linewidth=1.5)
legend("BPM", "LWPC", "FDTD");
ylabel("Amp (dB μV/m)");

dftdiff = namp - dft.amp
datdiff = namp .- dat.amp
ylim_min = floorto(min(quantile(dftdiff, 0.01), quantile(datdiff, 0.01)), 5) - 1
ylim_max = ceilto(max(quantile(dftdiff, 0.99), quantile(datdiff, 0.99)), 5) + 1

lowplot = subplot(4, 1, 4);
plot(dat.dist, datdiff, "r-", linewidth=1.5);
oplot(dat.dist, dftdiff, "y-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
ylabel("Δ BPM");
xlabel("Distance (km)");
yticks(5)
gcf()

savefig(joinpath(basepath,"homogeneous_amp.png"))

subplot(4, 1, (1, 3));
plot(X/1000, rad2deg.(phase), "b", linewidth=1.5)
oplot(X/1000, dat.phase, "r", linewidth=1.5)
oplot(X/1000, dft.phase, "y", linewidth=1.5)
legend("BPM", "LWPC", "FDTD");
ylabel("Phase (deg)");

dftdiff = rad2deg.(phase) - dft.phase
datdiff = rad2deg.(phase) .- dat.phase
ylim_min = floorto(min(quantile(dftdiff, 0.01), quantile(datdiff, 0.01)), 5) - 1
ylim_max = ceilto(max(quantile(dftdiff, 0.99), quantile(datdiff, 0.99)), 5) + 1

lowplot = subplot(4, 1, 4);
plot(dat.dist, datdiff, "r-", linewidth=1.5);
oplot(dat.dist, dftdiff, "y-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
ylabel("Δ BPM");
xlabel("Distance (km)");
yticks(fld(ylim_max-ylim_min,3))
gcf()

savefig(joinpath(basepath,"homogeneous_phase.png"))


# Homogeneous scenario 2
scenario = "homoscenario2"
raw = CSV.File(joinpath(basepath, scenario, scenario*".log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]

E, phase, amp = homoscenario2(3e9)
namp = copy(amp)

# TEMP  first point in amp is NaN which causes it to not plot
namp[1] = namp[2]

subplot(4, 1, (1, 3));
plot(X/1000, namp, linewidth=1.5)
oplot(X/1000, dat.amp, linewidth=1.5)
legend("BPM", "LWPC");
ylabel("Amp (dB μV/m)");

datdiff = namp-dat.amp
ylim_min = floor(quantile(datdiff, 0.01))
ylim_max = ceil(quantile(datdiff, 0.99))

lowplot = subplot(4, 1, 4);
plot(X/1000, datdiff, "r-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
yticks(1)
ylabel("Δ BPM");
xlabel("Distance (km)");
gcf()

savefig(joinpath(basepath,"homogeneous2.png"))




# Single transition
scenario = "singletransition"
raw = CSV.File(joinpath(basepath, scenario, "singletransition.log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]

rawdft = CSV.File(joinpath(basepath, scenario, "alldfts.csv");
                  header=false)
dft = DataFrame(dist=rawdft.Column1, amp=rawdft.Column2, phase=rawdft.Column3);
dft = dft[1:findfirst(dft.dist .== X[end]/1000),:]
dft = dft[1:10:end,:]
dft.amp .+= 116  # used for Jamesina comparison
dft.phase .-= 33

_, thresh_phase, thresh_amp = mc_scenario(3e9)
nthresh_amp = copy(thresh_amp)
nthresh_amp[1] = nthresh_amp[2]

# TODO: why?
thresh_phase .+= π

subplot(4, 1, (1, 3));
plot(X/1000, nthresh_amp, "b", linewidth=1.5);
oplot(X/1000, dat.amp, "r", linewidth=1.5);
oplot(X/1000, dft.amp, "y", linewidth=1.5);
legend("BPM_thresh", "LWPC", "FDTD");
ylabel("Amp (dB μV/m)");

dftdiff = nthresh_amp .- dft.amp
datdiff = nthresh_amp .- dat.amp
ylim_min = floorto(min(quantile(dftdiff, 0.01), quantile(datdiff, 0.01)), 1) - 1.1
ylim_max = ceilto(max(quantile(dftdiff, 0.99), quantile(datdiff, 0.99)), 1) + 1.1

lowplot = subplot(4, 1, 4);
plot(X/1000, datdiff, "r-", linewidth=1.5);
oplot(X/1000, dftdiff, "y-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
# ylim(-2, 2)
yticks(fld(ylim_max-ylim_min,3))
ylabel("Δ BPM");
xlabel("Distance (km)");
gcf()

savefig(joinpath(basepath,"single_amp.png"))


ylim_min = min(quantile(rad2deg.(thresh_phase), 0.01), quantile(dat.phase, 0.01), quantile(dft.phase, 0.01))
ylim_max = max(quantile(rad2deg.(thresh_phase), 0.99), quantile(dat.phase, 0.99), quantile(dft.phase, 0.99))

ylim_min = floorto(ylim_min, 5) - 1
ylim_max = ceilto(ylim_max, 5) + 1

subplot(4, 1, (1, 3));
plot(X/1000, rad2deg.(thresh_phase), "b", linewidth=1.5);
oplot(X/1000, dat.phase, "r", linewidth=1.5);
oplot(X/1000, dft.phase, "y", linewidth=1.5);
ylim(ylim_min, ylim_max)
legend("BPM_thresh", "LWPC", "FDTD", location="lower right");
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
yticks(fld(ylim_max-ylim_min,3))
xlabel("Distance (km)");
ylabel("Δ BPM");
gcf()

savefig(joinpath(basepath,"single_phase.png"))



# Double transition
scenario = "doubletransition"
raw = CSV.File(joinpath(basepath, scenario, "doubletransition.log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]

rawdft = CSV.File(joinpath(basepath, scenario, "alldfts.csv");
                  header=false)
dft = DataFrame(dist=rawdft.Column1, amp=rawdft.Column2, phase=rawdft.Column3);
dft = dft[1:findfirst(dft.dist .== X[end]/1000),:]
dft = dft[1:10:end,:]
dft.amp .+= 116  # used for Jamesina comparison
dft.phase .-= 33  # used for Jamesina comparison

_, thresh_phase, thresh_amp = mc_scenario_3(3e9)
nthresh_amp = copy(thresh_amp)
nthresh_amp[1] = nthresh_amp[2]

# TODO  why offset?
thresh_phase .+= π

subplot(4, 1, (1, 3));
plot(X/1000, nthresh_amp, "b", linewidth=1.5);
oplot(X/1000, dat.amp, "r", linewidth=1.5);
oplot(X/1000, dft.amp, "y", linewidth=1.5);
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
# ylim(-2, 2)
yticks(fld(ylim_max-ylim_min,3))
ylabel("Δ BPM");
xlabel("Distance (km)");
gcf()

savefig(joinpath(basepath,"double_amp.png"))


ylim_min = min(quantile(rad2deg.(thresh_phase), 0.01), quantile(dat.phase, 0.01), quantile(dft.phase, 0.01))
ylim_max = max(quantile(rad2deg.(thresh_phase), 0.99), quantile(dat.phase, 0.99), quantile(dft.phase, 0.99))

ylim_min = floorto(ylim_min, 5) - 1
ylim_max = ceilto(ylim_max, 5) + 1

# ylim_min = -130
# ylim_max = 380

subplot(4, 1, (1, 3));
plot(X/1000, rad2deg.(thresh_phase), "b", linewidth=1.5);
oplot(X/1000, dat.phase, "r", linewidth=1.5);
oplot(X/1000, dft.phase, "y", linewidth=1.5);
ylim(ylim_min, ylim_max)
legend("BPM_thresh", "LWPC", "FDTD", location="lower right");
ylabel("Phase (deg)");

dftdiff = rad2deg.(thresh_phase)-dft.phase
datdiff = rad2deg.(thresh_phase)-dat.phase
ylim_min = floorto(min(quantile(dftdiff, 0.01), quantile(datdiff, 0.01)), 5) - 1
ylim_max = ceilto(max(quantile(dftdiff, 0.99), quantile(datdiff, 0.99)), 5) + 1

# ylim_min = -20
# ylim_max = 20

subplot(4, 1, 4);
plot(dat.dist, datdiff, "r-", linewidth=1.5);
oplot(dat.dist, dftdiff, "y-", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylim(ylim_min, ylim_max)
yticks(ceilto(fld(ylim_max-ylim_min,3),2))
xlabel("Distance (km)");
ylabel("Δ BPM");
gcf()

savefig(joinpath(basepath,"double_phase.png"))
