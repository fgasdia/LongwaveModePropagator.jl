function test_modeterms(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    ea = TEST_MODES[scenario][1]
    modeequation = PhysicalModeEquation(ea, tx.frequency, waveguide)

    groundsampler = GroundSampler(rx.distance, rx.fieldcomponent)
    sampler = Sampler(rx.distance, rx.fieldcomponent, 0.0)

    # each field
    for fc in (Fields.Ez, Fields.Ey, Fields.Ex)
        groundsampler = GroundSampler(rx.distance, fc)
        sampler = Sampler(rx.distance, fc, 0.0)

        tx1, rx1 = LMP.modeterms(modeequation, tx, sampler)
        tx2, rx2 = LMP.modeterms(modeequation, tx, groundsampler)  # specialized

        @test tx1 ≈ tx2
        @test rx1 ≈ rx2
    end

    # frequency mismatch with modeequation
    txwrong = Transmitter(15e3)
    @test_throws ArgumentError LMP.modeterms(modeequation, txwrong, sampler)
end

function test_Efield(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    modes = TEST_MODES[scenario]

    X = LMP.distance(rx, tx)
    E = zeros(ComplexF64, length(X))

    E1 = LMP.Efield!(E, modes, waveguide, tx, rx)  # in-place
    E2 = LMP.Efield(modes, waveguide, tx, rx)  # out-of-place
    @test E1 == E2

    singlerx = GroundSampler(1000e3, rx.fieldcomponent)
    E3 = LMP.Efield(modes, waveguide, tx, singlerx)  # specialized
    distidx = findfirst(x->x==1000e3, X)
    @test E3 ≈ E2[distidx]

    # E and X length mismatch
    Ewrong = zeros(ComplexF64, 13)
    @test_throws ArgumentError LMP.Efield!(Ewrong, modes, waveguide, tx, rx)
end


@testset "modesum.jl" begin
    @info "Testing modesum"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        test_modeterms(scn)
        test_Efield(scn)
    end
end
