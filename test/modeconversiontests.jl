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


struct Wavefields{T,T2}
    # fields(z) for each mode, v[mode][height]
    v::Vector{Vector{SVector{6,T}}}
    eas::Vector{EigenAngle}
    zs::T2
end

function Wavefields(eas::Vector{EigenAngle}, zs::T2) where {T2<:AbstractVector}
    Wavefields{T,T2}([Vector{SVector{6,T}}(undef,length(zs)) for i = 1:length(eas)], eas, zs)
end

Base.getindex(A::Wavefields, i::Int) = A.v[i]

eigenangles(A::Wavefields) = A.eas
heights(A::Wavefields) = A.zs
numeigenangles(A::Wavefields) = length(A.eas)
numheights(A::Wavefields) = length(A.zs)



function test_calculate_wavefields!()
    bfield, tx, ground, electrons, ea, zs = resonant_scenario()

    wavefields = Wavefields(ea, zs)
    adjoint_wavefields = Wavefields(ea, zs)

    calculate_wavefields!(wavefields, adjoint_wavefields,
                          bfield, tx.frequency, ground, electrons)

    return wavefields, adjoint_wavefields
end


function mc_scenario()
    waveguide = HomogeneousWaveguide[]

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(68), deg2rad(111)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 75, 0.32),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001)))

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, deg2rad(68), deg2rad(111)),
                                          Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 78, 0.35),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001)))

    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
    rx = GroundSampler(0:10e3:2000e3, LWMS.FC_Ez)



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



function lwpce(dst, waveguide, tx, rx)
    frequency = tx.frequency

    k = frequency.k

    mik = complex(0, -k)
    aconst = -8686*k
    econst = 20log10(35*k)
    sum0 = 682.2408*sqrt(frequency.f*tx.power)

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


    nsgmnt = 1
    nrsgmnt = length(waveguide)

    modes = LWMS.findmodes(origcoords, frequency, waveguide[nsgmnt], tolerance)

    wavefields = Wavefields(modes, zs)
    adjwavefields = Wavefields(modes, zs)

    calculate_wavefields!(wavefields, adjwavefields, frequency, waveguide[nsgmnt])

    # soln_a is for `Hy`
    soln_b = Vector{ComplexF64}(undef, nreigen2)

    dst0 = 99999
    if dst < dst0
        xone = 0
        if nrgsgmnt == 1
            xtwo = 40000
        else
            xtwo = waveguide[nsgmnt+1].distance
        end

        wvg = waveguide[nsgmnt]
        eas = eigenangles(wvg)
        nreigen2 = length(eas)

        # `for` loop from 180 to 193 and LW_STEP and LW_HTGAIN all contained in modeterms
        for m2 = 1:nreigen2
            ta, _ = modeterms(eas[m2], frequency, wvg, emitter_orientation, sampler_orientation)
            soln_b[m2] = ta
        end
    end

    temp = Vector{ComplexF64}(undef, nreigen2)
    while dst > xtwo
        # End of current slab
        mikx = mik*(xtwo - xone)
        nreigen1 = nreigen2
        for m1 = 1:nreigen1
            S₀ = referencetoground(eas[m1].sinθ)
            # Exctation factors at end of slab. LWPC uses `Hy`
            temp[m1] = soln_b[m1]*exp(mikx*(S₀ - 1))
        end
        xone = xtwo

        # Load next slab
        nsgmnt += 1
        wvg = waveguide[nsgmnt]

        if nsgmnt < nrsgmnt
            xtwo = waveguide[nsgmnt+1].distance
        else
            xtwo = 40000
        end

        a = modeconversion(prevwavefields, wavefields, adjwavefields)
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
        S₀ = referencetoground(eas[m2].sinθ)
        tb += soln_b[m2]*exp(mikx*(S₀ - 1))*factor
    end

    dst0 = dst

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
