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

    modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8

    modes = LWMS.findmodes(origcoords, modeparams, tolerance)
    # ea = modes[argmax(real(modes))]  # largest real resonant mode

    ea = [LWMS.EigenAngle(m) for m in modes]

    return bfield, tx, ground, electrons, ea, zs
end


struct Wavefields{T,T2<:AbstractVector}
    # fields(z) for each mode, v[mode][height]
    # TODO: Are wavefields always complex (independent of eas?)
    v::Vector{Vector{SVector{6,T}}}
    eas::Vector{LWMS.EigenAngle{T}}
    zs::T2
end

function Wavefields(eas::Vector{LWMS.EigenAngle{T}},zs::T2) where {T, T2<:AbstractVector}
    Wavefields{T,T2}([Vector{SVector{6,T}}(undef,length(zs)) for i = 1:length(eas)], eas, zs)
end

Base.getindex(A::Wavefields, i::Int) = A.v[i]

eigenangles(A::Wavefields) = A.eas
heights(A::Wavefields) = A.zs
numeigenangles(A::Wavefields) = length(A.eas)
numheights(A::Wavefields) = length(A.zs)

function calculate_wavefields!(wavefields, adjoint_wavefields,
                               bfield, frequency, ground, species) where {T}

    @assert heights(wavefields) == heights(adjoint_wavefields)

    zs = heights(wavefields)
    modes = eigenangles(wavefields)

    @inbounds for m in eachindex(modes)
        EH = LWMS.fieldstrengths(zs, modes[m], frequency, bfield, species, ground)

        adjoint_bfield = BField(bfield.B, -bfield.dcl, bfield.dcm, bfield.dcn)
        EHadjoint = LWMS.fieldstrengths(zs, modes[m], frequency,adjoint_bfield, species, ground)

        copyto!(wavefields[m], EH)
        copyto!(adjoint_wavefields[m], EHadjoint)
    end

    return nothing
end

function test_calculate_wavefields!()
    bfield, tx, ground, electrons, ea, zs = resonant_scenario()

    wavefields = Wavefields(ea, zs)
    adjoint_wavefields = Wavefields(ea, zs)

    calculate_wavefields!(wavefields, adjoint_wavefields,
                          bfield, tx.frequency, ground, electrons)

    return wavefields, adjoint_wavefields
end



function mc_scenario()
    waveguide = LWMS.WaveguideSegment[]

    push!(waveguide, LWMS.WaveguideSegment(BField(50e-6, deg2rad(68), deg2rad(111)),
                                      Species(qₑ, mₑ,
                                                  z -> waitprofile(z, 75, 0.32),
                                                  electroncollisionfrequency),
                                      Ground(15, 0.001)))

    push!(waveguide, LWMS.WaveguideSegment(BField(50e-6, deg2rad(68), deg2rad(111)),
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
        modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, species)

        origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
        origcoords .= deg2rad.(origcoords)
        tolerance = 1e-8

        modes = LWMS.findmodes(origcoords, modeparams, tolerance)
        ea = [LWMS.EigenAngle(m) for m in modes]

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

"""
this function takes just 1 mode conversion step
"""
function modeconversion(previous_wavefields::Wavefields{T,T2},
                        wavefields::Wavefields{T,T2}, adjwavefields::Wavefields{T,T2},
                        transmitter_slab=false) where {T,T2}

    @assert numheights(previous_wavefields) == numheights(wavefields)
    zs = heights(wavefields)

    # TODO: Assuming `length(zs)` is always the same, we can reuse `product`
    product = Vector{T}(undef, numheights(wavefields))

    # Calculate normalization terms
    N = Vector{T}(undef, numeigenangles(wavefields))
    for m in eachindex(eigenangles(wavefields))
        for i in eachindex(zs)
            # TODO: is a view faster?
            @inbounds f = wavefields[m][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
            @inbounds g = SVector{4,T}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
            product[i] = adjoint(g)*f
        end

        N[m] = integrate(zs, product, RombergEven())
        @test isapprox(N[m], trapz(zs, product), rtol=1e-3)
        @test isapprox(N[m], romberg(zs, product), rtol=1e-3)
    end

    if transmitter_slab
        a = zeros(T, length(eas), length(eas))
        for i in diagind(a)
            @inbounds a[i] = one(ComplexF64)
        end

        # TODO: a = I ?  # UniformScaling identity matrix, probably need to pull this out for type consistency
    else
        I = Matrix{T}(undef, numeigenangles(adjwavefields), numeigenangles(previous_wavefields))
        for k in eachindex(eigenangles(previous_wavefields))
            for m in eachindex(eigenangles(adjwavefields))
                for i in eachindex(heights(adjwavefields))
                    @inbounds f = previous_wavefields[k][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
                    @inbounds g = SVector{4,T}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
                    product[i] = adjoint(g)*f
                end
                I[m,k] = integrate(zs, product, RombergEven())
            end
        end

        # Total conversion
        for k in eachindex(eigenangles(previous_wavefields))
            for m in eachindex(eigenangles(adjwavefields))
                a[k,m] = I[m,k]/N[m]
            end
        end
    end

    return a
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
