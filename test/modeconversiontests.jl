using Test
using LinearAlgebra
using StaticArrays
using Plots
using NumericalIntegration
using Trapz  # for testing only
using Parameters

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
    waveguide = HomogeneousWaveguide[]

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

    waveguide = LWMS.SegmentedWaveguide(waveguide)

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


    # ztop = LWMS.TOPHEIGHT
    # zs = range(ztop, zero(ztop), length=257)
    #
    # waveguide_wavefields = Wavefields[]
    # waveguide_adjwavefields = similar(waveguide_wavefields)
    # for s in eachindex(waveguide)
    #     @unpack bfield, species, ground = waveguide[s]
    #
    #     origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    #     origcoords .= deg2rad.(origcoords)
    #     tolerance = 1e-8
    #
    #     ea = LWMS.findmodes(origcoords, tx.frequency, waveguide[s], tolerance)
    #
    #     wavefields = Wavefields(ea, zs)
    #     adjwavefields = Wavefields(ea, zs)
    #
    #     calculate_wavefields!(wavefields, adjwavefields,
    #                           bfield, tx.frequency, ground, species)
    #
    #     # TODO: only store previous and make sure size reflects length(ea)
    #     push!(waveguide_wavefields, wavefields)
    #     push!(waveguide_adjwavefields, adjwavefields)
    # end
    #
    # # Try mode conversion
    # modeconversion(waveguide_wavefields[1],
    #                waveguide_wavefields[2], waveguide_adjwavefields[2])
end


basepath = "/home/forrest/research/LAIR/ModeSolver/lwpc_comparisons/"

raw = CSV.File(joinpath(basepath, "singletransition.log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]

widedf = DataFrame(dist=dat.dist,
                   lwpc_amp=dat.amp, lwpc_phase=dat.phase,
                   lwms_amp=amp, lwms_phase=rad2deg.(phase))

CSV.write(joinpath(basepath, "singletransition.csv"), widedf)

outpath = "/home/forrest/UCB/SP_2020/PropagationModeling/figures"

function myplot(bfield, ground, tx)
    f = tx.frequency.f/1000
    B = Int(bfield.B/1e-9)
    dip = Int(rad2deg(LWMS.dip(bfield)))
    az = Int(rad2deg(LWMS.azimuth(bfield)))
    epsr = ground.ϵᵣ
    sigma = ground.σ

    gp_title = """
    TITLE = '"$f kHz\\n\\
    |B|: $B nT, dip: $(dip)°, az: $(az)°\\n\\
    h\'\': 75 - 70, β: 0.32 - 0.25 \\n\\
    ϵ_r: $epsr, σ: $sigma"'"""

    open("gp_title", "w") do io
        write(io, gp_title)
    end

    ga = `gnuplot -c "$(joinpath(outpath,"vertical_amp_linux.gp"))" "$(joinpath(basepath,"singletransition.csv"))" "$(joinpath(outpath,""))"`
    run(ga)

    gp = `gnuplot -c "$(joinpath(outpath,"vertical_phase_linux.gp"))" "$(joinpath(basepath,"singletransition.csv"))" "$(joinpath(outpath,""))"`
    run(gp)
end
bfield = waveguide[1].bfield
ground = waveguide[1].ground
myplot(bfield, ground, tx)






function lwpce(waveguide, tx, rx)
    X = LWMS.distance(rx, tx)
    E = Vector{ComplexF64}(undef, length(X))

    frequency = tx.frequency

    k = frequency.k

    mik = complex(0, -k)
    sum0 = 682.2408*sqrt(frequency.f/1000*tx.power/1000)

    # Antenna orientation factors
    Sγ, Cγ = sincos(π/2 - LWMS.elevation(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(LWMS.azimuth(tx))  # ϕ is measured from `x`
    t1, t2, t3 = Cγ, Sγ*Sϕ, Sγ*Cϕ

    zt = LWMS.altitude(tx)

    rxcomponent = LWMS.fieldcomponent(rx)
    zr = LWMS.altitude(rx)

    emitter_orientation = (Sγ=Sγ, Cγ=Cγ, Sϕ=Sϕ, Cϕ=Cϕ, zt=zt)
    sampler_orientation = (rxcomponent=rxcomponent, zr=zr)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8
    zs = range(LWMS.TOPHEIGHT, 0.0, length=257)


    nrsgmnt = length(waveguide)

    wavefields_vec = Vector{LWMS.Wavefields}(undef, nrsgmnt)
    adjwavefields_vec = Vector{LWMS.Wavefields}(undef, nrsgmnt)

    # TEMP Precalculate all wavefields
    for nsgmnt in 1:nrsgmnt
        modes = LWMS.findmodes(origcoords, frequency, waveguide[nsgmnt], tolerance)

        wavefields = LWMS.Wavefields(modes, zs)
        adjwavefields = LWMS.Wavefields(modes, zs)

        LWMS.calculate_wavefields!(wavefields, adjwavefields, frequency, waveguide[nsgmnt])

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

        xone = zero(wvg.distance)
        # if nrsgmnt == 1
        #     xtwo = 40000e3
        # else
        xtwo = waveguide[nsgmnt+1].distance

        # soln_a is for `Hy`
        soln_b = Vector{ComplexF64}(undef, nreigen2)

        # `do` loop from 180 to 193 and LW_STEP and LW_HTGAIN all contained in modeterms
        for m2 = 1:nreigen2
            ta, tb = LWMS.modeterms(eas[m2], frequency, wvg, emitter_orientation, sampler_orientation)
            soln_b[m2] = ta*tb
        end


        temp = Vector{ComplexF64}(undef, nreigen2)
        while dst > xtwo
            # End of current slab
            mikx = mik*(xtwo - xone)
            nreigen1 = nreigen2
            for m1 = 1:nreigen1
                S₀ = LWMS.referencetoground(eas[m1].sinθ)
                # Exctation factors at end of slab. LWPC uses `Hy`
                temp[m1] = soln_b[m1]*exp(mikx*(S₀ - 1))
            end
            xone = xtwo

            # Load next slab
            nsgmnt += 1
            wvg = waveguide[nsgmnt]

            wavefields = wavefields_vec[nsgmnt]
            adjwavefields = adjwavefields_vec[nsgmnt]
            prevwavefields = wavefields_vec[nsgmnt-1]

            nreigen2 = LWMS.numeigenangles(wavefields)

            if nsgmnt < nrsgmnt
                xtwo = waveguide[nsgmnt+1].distance
            else
                xtwo = oftype(wvg.distance, 40_000_000)
            end

            a = LWMS.modeconversion(prevwavefields, wavefields, adjwavefields)

            soln_b = zeros(ComplexF64, nreigen2)
            for m2 = 1:nreigen2
                for m1 = 1:nreigen1
                    soln_b[m2] += temp[m1]*a[m1,m2]
                end
                # ta is `Hy` excitation factor
                # Then LWPC calculates E into soln_b
            end

        end

        mikx = mik*(dst - xone)
        factor = sum0/sqrt(abs(sin(dst/EARTH_RADIUS)))

        # For each component
        tb = zero(ComplexF64)
        for m2 = 1:nreigen2
            S₀ = LWMS.referencetoground(eas[m2].sinθ)
            tb += soln_b[m2]*exp(mikx*(S₀ - 1))*factor
        end

        E[i] = tb
    end

    return E
end













function lwpce(waveguide, tx, rx)
    X = LWMS.distance(rx, tx)
    E = Vector{ComplexF64}(undef, length(X))

    frequency = tx.frequency

    k = frequency.k

    mik = complex(0, -k)
    sum0 = 682.2408*sqrt(frequency.f/1000*tx.power/1000)

    # Antenna orientation factors
    Sγ, Cγ = sincos(π/2 - LWMS.elevation(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(LWMS.azimuth(tx))  # ϕ is measured from `x`
    t1, t2, t3 = Cγ, Sγ*Sϕ, Sγ*Cϕ

    zt = LWMS.altitude(tx)

    rxcomponent = LWMS.fieldcomponent(rx)
    zr = LWMS.altitude(rx)

    emitter_orientation = (Sγ=Sγ, Cγ=Cγ, Sϕ=Sϕ, Cϕ=Cϕ, zt=zt)
    sampler_orientation = (rxcomponent=rxcomponent, zr=zr)

    # const = sum0

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8
    zs = range(LWMS.TOPHEIGHT, 0.0, length=257)

    nrsgmnt = length(waveguide)
    for i in eachindex(X)
        dst = X[i]

        nsgmnt = 1
        modes = LWMS.findmodes(origcoords, frequency, waveguide[nsgmnt], tolerance)

        wavefields = LWMS.Wavefields(modes, zs)
        adjwavefields = LWMS.Wavefields(modes, zs)
        prevwavefields = LWMS.Wavefields(modes, zs)

        LWMS.calculate_wavefields!(wavefields, adjwavefields, frequency, waveguide[nsgmnt])

        nreigen2 = LWMS.numeigenangles(wavefields)

        # soln_a is for `Hy`
        soln_b = Vector{ComplexF64}(undef, nreigen2)

        # dst0 = 99999 the dst vs dst0 checks aren't really needed because dist issorted
        xone = 0
        if nrsgmnt == 1
            xtwo = 40000e3
        else
            xtwo = waveguide[nsgmnt+1].distance
        end

        wvg = waveguide[nsgmnt]
        eas = LWMS.eigenangles(wavefields)

        # `for` loop from 180 to 193 and LW_STEP and LW_HTGAIN all contained in modeterms
        for m2 = 1:nreigen2
            ta, _ = LWMS.modeterms(eas[m2], frequency, wvg, emitter_orientation, sampler_orientation)
            soln_b[m2] = ta
        end


        temp = Vector{ComplexF64}(undef, nreigen2)
        while dst > xtwo
            # End of current slab
            mikx = mik*(xtwo - xone)
            nreigen1 = nreigen2
            for m1 = 1:nreigen1
                S₀ = LWMS.referencetoground(eas[m1].sinθ)
                # Exctation factors at end of slab. LWPC uses `Hy`
                temp[m1] = soln_b[m1]*exp(mikx*(S₀ - 1))
            end
            xone = xtwo

            # Load next slab
            prevwavefields, wavefields = wavefields, prevwavefields
            nsgmnt += 1
            wvg = waveguide[nsgmnt]

            LWMS.calculate_wavefields!(wavefields, adjwavefields, frequency, wvg)

            if nsgmnt < nrsgmnt
                xtwo = waveguide[nsgmnt+1].distance
            else
                xtwo = 40000e3
            end

            a = LWMS.modeconversion(prevwavefields, wavefields, adjwavefields)

            soln_b = zeros(ComplexF64, nreigen2)
            for m2 = 1:nreigen2
                for m1 = 1:nreigen1
                    soln_b[m2] += temp[m1]*a[m1,m2]
                end
                # ta is `Hy` excitation factor
                # Then LWPC calculates E into soln_b
            end

        end

        mikx = mik*(dst - xone)
        factor = sum0/sqrt(abs(sin(dst/EARTH_RADIUS)))

        # For each component
        tb = zero(ComplexF64)
        for m2 = 1:nreigen2
            S₀ = LWMS.referencetoground(eas[m2].sinθ)
            tb += soln_b[m2]*exp(mikx*(S₀ - 1))*factor
        end

        E[i] = tb
    end

    return E
end










i = 6
EH = reshape(reinterpret(ComplexF64, EH), 6, :)
plot(real(EH[i,:]), zs/1000, color="black")
plot!(imag(EH[i,:]), zs/1000, color="black")
plot!(abs.(EH[i,:]), zs/1000, color="black")
plot!(-abs.(EH[i,:]), zs/1000, color="black")

EHadjoint = reshape(reinterpret(ComplexF64, EHadjoint), 6, :)
plot!(real(EHadjoint[i,:]), zs/1000, color="blue")
plot!(imag(EHadjoint[i,:]), zs/1000, color="blue")
plot!(abs.(EHadjoint[i,:]), zs/1000, color="blue")
plot!(-abs.(EHadjoint[i,:]), zs/1000, color="blue")


x = 0.28656210625768974
@test LWMS.pow23(x) == LWMS.pow23(complex(x))
y = 100*x
@test LWMS.pow23(y) == LWMS.pow23(complex(y))





using GRUtils

basepath = "/home/forrest/research/LAIR/ModeSolver/lwpc_comparisons/"

raw = CSV.File(joinpath(basepath, "singletransition.log");
               skipto=40, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

dat = dat[1:401,:]


_, phase, amp = mc_scenario()

X = LWMS.distance(rx,tx);
namp = replace(amp, NaN=>0.0)
plot(X/1000, namp, linewidth=2)
oplot(dat.dist, dat.amp)
xlabel("Distance (km)")
ylabel("Amp")

plot(X/1000, rad2deg.(phase), linewidth=2)
oplot(dat.dist, dat.phase)
xlabel("Distance (km)")
ylabel("Phase")
