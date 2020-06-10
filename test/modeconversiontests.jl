using Test
using LinearAlgebra
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


function resonant_scenario()
    bfield = BField(50e-6, deg2rad(68), deg2rad(111))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
    ground = Ground(15, 0.001)
    electrons = Species(qₑ, mₑ,
                            z -> waitprofile(z, 75, 0.32),
                            electroncollisionfrequency)

    ztop = LWMS.TOPHEIGHT
    zs = range(ztop, zero(ztop), length=257)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8

    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    modes = LWMS.findmodes(origcoords, tx.frequency, waveguide, tolerance)
    # ea = modes[argmax(real(modes))]  # largest real resonant mode

    return bfield, tx, ground, electrons, ea, zs
end



function test_calculate_wavefields!()
    bfield, tx, ground, electrons, ea, zs = resonant_scenario()

    wavefields = Wavefields(ea, zs)
    adjoint_wavefields = Wavefields(ea, zs)

    calculate_wavefields!(wavefields, adjoint_wavefields,
                          bfield, tx.frequency, ground, electrons)

    return wavefields, adjoint_wavefields
end


function homoscenario()
    waveguide =  HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 75, 0.32),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001))

    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    rx = GroundSampler(0:5e3:2000e3, LWMS.FC_Ez)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8

    modes = LWMS.findmodes(origcoords, tx.frequency, waveguide, tolerance)

    LWMS.Efield(modes, waveguide, tx, rx)
end


function mc_scenario()
    waveguide = LWMS.SegmentedWaveguide(HomogeneousWaveguide)

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 75, 0.32),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001)))

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(90), deg2rad(0)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 70, 0.25),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001), 1000e3))


    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
    rx = GroundSampler(0:5e3:2000e3, LWMS.FC_Ez)


    E = lwpce(waveguide, tx, rx)

    phase = Vector{Float64}(undef, length(E))
    amp = Vector{Float64}(undef, length(E))
    @inbounds for i in eachindex(E)
        e = E[i]
        phase[i] = angle(e)  # ranges between -π:π rad
        amp[i] = 10log10(abs2(e))  # == 20log10(abs(E))
    end

    # By definition, phase at transmitter is 0, but is calculated as NaN
    if isnan(phase[1])
        phase[1] = 0
    end

    LWMS.unwrap!(phase)

    return E, phase, amp
end


function lwpce(waveguide, tx, rx)
    X = LWMS.distance(rx, tx)
    E = Vector{ComplexF64}(undef, length(X))

    frequency = tx.frequency
    k = frequency.k

    sum0 = 682.2408*sqrt(frequency.f/1000*tx.power/1000)

    # Antenna orientation factors
    Sγ, Cγ = sincos(π/2 - LWMS.elevation(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(LWMS.azimuth(tx))  # ϕ is measured from `x`

    zt = LWMS.altitude(tx)

    rxcomponent = LWMS.fieldcomponent(rx)
    zr = LWMS.altitude(rx)

    emitter_orientation = (t1=Cγ, t2=Sγ*Sϕ, t3=Sγ*Cϕ, zt=zt)
    sampler_orientation = (rxcomponent=rxcomponent, zr=zr)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8
    zs = range(LWMS.TOPHEIGHT, 0.0, length=513)


    nrsgmnt = length(waveguide)

    wavefields_vec = Vector{LWMS.Wavefields{typeof(zs)}}(undef, nrsgmnt)
    adjwavefields_vec = Vector{LWMS.Wavefields{typeof(zs)}}(undef, nrsgmnt)

    # TEMP Precalculate all wavefields
    for nsgmnt in 1:nrsgmnt
        wvg = waveguide[nsgmnt]
        modes = LWMS.findmodes(origcoords, frequency, wvg, tolerance)

        # adjoint wavefields are wavefields through adjoint waveguide, but for same modes as wavefield
        @unpack bfield, species, ground = wvg
        adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
        adjwvg = LWMS.HomogeneousWaveguide(adjoint_bfield, species, ground)

        wavefields = LWMS.Wavefields(modes, zs)
        adjwavefields = LWMS.Wavefields(modes, zs)

        LWMS.calculate_wavefields!(wavefields, adjwavefields, frequency, wvg, adjwvg)

        wavefields_vec[nsgmnt] = wavefields
        adjwavefields_vec[nsgmnt] = adjwavefields
    end

    for i in eachindex(X)
        dst = X[i]

        nsgmnt = 1

        wvg = waveguide[nsgmnt]
        wavefields = wavefields_vec[nsgmnt]
        eas = LWMS.eigenangles(wavefields)
        nreigen2 = length(eas)

        xone = wvg.distance  # @assert wvg.distance == 0
        xtwo = waveguide[nsgmnt+1].distance

        # soln_a is for `Hy`
        soln_b = Vector{ComplexF64}(undef, nreigen2)

        # `do` loop from 180 to 193 and LW_STEP and LW_HTGAIN all contained in modeterms
        for m2 = 1:nreigen2
            ta, tb = LWMS.modeterms(eas[m2], frequency, wvg, emitter_orientation, sampler_orientation)
            soln_b[m2] = ta*tb
        end

        while dst > xtwo
            # End of current slab
            x = xtwo - xone
            mkx = -k*x
            temp = Vector{ComplexF64}(undef, nreigen2)
            for m2 = 1:nreigen2
                S₀ = LWMS.referencetoground(eas[m2].sinθ)
                # Excitation factors at end of slab. LWPC uses `Hy`
                soln_b[m2] *= cis(mkx*(S₀ - 1))
                temp[m2] = soln_b[m2]
            end
            nreigen1 = nreigen2
            xone = xtwo

            # Load next slab
            nsgmnt += 1
            wvg = waveguide[nsgmnt]

            wavefields = wavefields_vec[nsgmnt]
            adjwavefields = adjwavefields_vec[nsgmnt]
            prevwavefields = wavefields_vec[nsgmnt-1]
            eas = LWMS.eigenangles(wavefields)
            nreigen2 = length(eas)

            if nsgmnt < nrsgmnt
                xtwo = waveguide[nsgmnt+1].distance
            else
                xtwo = oftype(wvg.distance, 40_000_000)
            end

            a = LWMS.modeconversion(prevwavefields, wavefields, adjwavefields)

            soln_b = a*temp
        end

        x = dst - xone
        mkx = -k*x
        factor = sum0/sqrt(abs(sin(dst/EARTH_RADIUS)))

        tb = zero(ComplexF64)
        for m2 = 1:nreigen2
            S₀ = LWMS.referencetoground(eas[m2].sinθ)
            tb += soln_b[m2]*cis(mkx*(S₀ - 1))
        end

        E[i] = tb*factor
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
colorscheme("dark")  # matches atom, otherwise "light"

basepath = "/home/forrest/research/LAIR/ModeSolver/lwpc_comparisons/"

raw = CSV.File(joinpath(basepath, "singletransition.log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]


_, phase, amp = mc_scenario()

X = LWMS.distance(rx,tx);
namp = copy(amp)

# TEMP  first point in amp is NaN which causes it to not plot
namp[1] = dat.amp[1]

subplot(4, 1, (1, 3));
plot(X/1000, namp, linewidth=1.5);
oplot(dat.dist, dat.amp, linewidth=1.5);
legend("BPM", "LWPC");
ylabel("Amp (dB μV/m)");

lowplot = subplot(4, 1, 4);
plot(dat.dist, namp-dat.amp, "gray", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
ylabel("BPM − LWPC");
xlabel("Distance (km)");
GRUtils.gcf()


subplot(4, 1, (1, 3));
plot(X/1000, rad2deg.(phase), linewidth=1.5);
oplot(dat.dist, dat.phase, linewidth=1.5);
legend("BPM", "LWPC");
ylabel("Phase (deg)");

subplot(4, 1, 4);
plot(dat.dist, rad2deg.(phase)-dat.phase, "gray", linewidth=1.5);
lowplot.viewport.inner[3] = 0.075;  # override bottom margin
lowplot.viewport.inner[4] = 0.2;  # and top margin
xlabel("Distance (km)");
ylabel("BPM − LWPC");
GRUtils.gcf()
