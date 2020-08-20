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


function bpm(waveguide::SegmentedWaveguide, tx, rx)
    zs = range(LWMS.TOPHEIGHT, 0, length=513)
    # zs = range(TOPHEIGHT, 0, length=129)
    nrsgmnt = length(waveguide)

    wavefields_vec = Vector{LWMS.Wavefields{typeof(zs)}}(undef, nrsgmnt)
    adjwavefields_vec = Vector{LWMS.Wavefields{typeof(zs)}}(undef, nrsgmnt)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-7

    for nsgmnt in 1:nrsgmnt
        wvg = waveguide[nsgmnt]
        modes = LWMS.findmodes(origcoords, tx.frequency, wvg, tolerance)

        # adjoint wavefields are wavefields through adjoint waveguide, but for same modes as wavefield
        @unpack bfield, species, ground = wvg
        adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
        adjwvg = HomogeneousWaveguide(adjoint_bfield, species, ground)

        wavefields = LWMS.Wavefields(modes, zs)
        adjwavefields = LWMS.Wavefields(modes, zs)

        LWMS.calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg)

        wavefields_vec[nsgmnt] = wavefields
        adjwavefields_vec[nsgmnt] = adjwavefields
    end

    E = Efield(waveguide, wavefields_vec, adjwavefields_vec, tx, rx)

    phase = Vector{Float64}(undef, length(E))
    amp = similar(phase)
    @inbounds for i in eachindex(E)
        e = E[i]
        phase[i] = angle(e)  # ranges between -π:π rad
        amp[i] = 10log10(abs2(e))  # == 20log10(abs(E))
    end

    # By definition, phase at transmitter is 0, but is calculated as NaN
    isnan(phase[1]) && (phase[1] = 0)
    unwrap!(phase)

    return E, phase, amp
end


function Efield(waveguide, wavefields_vec, adjwavefields_vec, tx, rx)
    X = distance(rx, tx)
    maxX = maximum(X)
    Xlength = length(X)
    E = Vector{ComplexF64}(undef, Xlength)

    frequency = tx.frequency
    k = frequency.k

    sum0 = 682.2408*sqrt(frequency.f/1000*tx.power/1000)

    # Antenna orientation factors
    Sγ, Cγ = sincos(π/2 - elevation(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(azimuth(tx))  # ϕ is measured from `x`

    zt = altitude(tx)

    rxcomponent = fieldcomponent(rx)
    zr = altitude(rx)

    emitter_orientation = (t1=Cγ, t2=Sγ*Sϕ, t3=Sγ*Cϕ, zt=zt)
    sampler_orientation = (rxcomponent=rxcomponent, zr=zr)

    nrsgmnt = length(waveguide)

    # initialize
    nreigen1 = 0
    temp = Vector{ComplexF64}(undef, 0)
    soln_a = similar(temp)
    soln_b = similar(temp)

    i = 1
    for nsgmnt = 1:nrsgmnt
        wvg = waveguide[nsgmnt]
        wavefields = wavefields_vec[nsgmnt]
        eas = eigenangles(wavefields)
        nreigen2 = length(eas)

        xone = wvg.distance  # @assert wvg.distance == 0 for nsgmnt == 1
        maxX < xone && break  # waveguide extends beyond greatest distance

        if nsgmnt < nrsgmnt
            xtwo = waveguide[nsgmnt+1].distance
        else
            xtwo = typemax(typeof(xone))
        end

        # soln_a is for `Hy`
        resize!(soln_a, nreigen2)
        resize!(soln_b, nreigen2)
        if nsgmnt > 1
            adjwavefields = adjwavefields_vec[nsgmnt]
            prevwavefields = wavefields_vec[nsgmnt-1]
            a = modeconversion(prevwavefields, wavefields, adjwavefields)
        end

        for m2 = 1:nreigen2
            ta, tb = modeterms(eas[m2], frequency, wvg, emitter_orientation, sampler_orientation)
            if nsgmnt == 1
                # Transmitter height gain only needed in transmitter slab
                soln_a[m2] = ta
            else
                # Otherwise, mode conversion
                soln_a_sum = zero(eltype(soln_a))
                for m1 = 1:nreigen1
                    soln_a_sum += temp[m1]*a[m2,m1]  # TODO: swap a m2,m1?
                end
                soln_a[m2] = soln_a_sum
            end

            soln_b[m2] = soln_a[m2]*tb
        end

        while X[i] < xtwo
            x = X[i] - xone
            factor = sum0/sqrt(abs(sin(X[i]/EARTH_RADIUS)))

            tb = zero(ComplexF64)
            for m2 = 1:nreigen2
                S₀ = referencetoground(eas[m2].sinθ)
                # tb += soln_b[m2]*cis(-k*x*(S₀))
                tb += soln_b[m2]*cis(-k*x*(S₀ - 1))
            end

            E[i] = tb*factor
            i += 1
            i > Xlength && break
        end

        if nsgmnt < nrsgmnt
            # End of current slab
            x = xtwo - xone

            resize!(temp, nreigen2)
            for m2 = 1:nreigen2
                S₀ = referencetoground(eas[m2].sinθ)
                # Excitation factors at end of slab. LWPC uses `Hy`
                soln_a[m2] *= cis(-k*x*(S₀ - 1))
                # soln_a[m2] *= cis(-k*x*(S₀))
                temp[m2] = soln_a[m2]
            end
            nreigen1 = nreigen2
        end
    end

    return E
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

savefig("homogeneous_amp.png")

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

savefig("homogeneous_phase.png")


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

savefig("homogeneous2.png")




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

savefig("single_amp.png")


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
xlabel("Distance (km)");
ylabel("Δ BPM");
gcf()

savefig("single_phase.png")



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

savefig("double_amp.png")


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

savefig("double_phase.png")
i
