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


struct Wavefields{T,T2<:AbstractVector}
    # fields(z) for each mode, v[mode][height]
    v::Vector{Vector{SVector{6,complex(T)}}}
    eas::Vector{EigenAngle{T}}
    zs::T2
end

function Wavefields(eas::Vector{EigenAngle{T}}, zs::T2) where {T,T2<:AbstractVector}
    Wavefields{T,T2}([Vector{SVector{6,complex(T)}}(undef,length(zs)) for i = 1:length(eas)], eas, zs)
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

    ztop = LWMS.TOPHEIGHT
    zs = range(ztop, zero(ztop), length=257)

    waveguide_wavefields = Wavefields[]
    waveguide_adjwavefields = similar(waveguide_wavefields)
    for s in eachindex(waveguide)
        @unpack bfield, species, ground = waveguide[s]

        origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
        origcoords .= deg2rad.(origcoords)
        tolerance = 1e-8

        ea = LWMS.findmodes(origcoords, tx.frequency, waveguide[s], tolerance)

        wavefields = Wavefields(ea, zs)
        adjwavefields = Wavefields(ea, zs)

        calculate_wavefields!(wavefields, adjwavefields,
                              bfield, tx.frequency, ground, species)

        # TODO: only store previous and make sure size reflects length(ea)
        push!(waveguide_wavefields, wavefields)
        push!(waveguide_adjwavefields, adjwavefields)
    end

    # Try mode conversion
    modeconversion(waveguide_wavefields[1],
                   waveguide_wavefields[2], waveguide_adjwavefields[2])
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
