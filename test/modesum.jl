function test_modeterms(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    ea = TEST_MODES[scenario][1]
    modeequation = PhysicalModeEquation(ea, tx.frequency, waveguide)

    groundsampler = GroundSampler(rx.distance, rx.fieldcomponent)
    sampler = Sampler(rx.distance, rx.fieldcomponent, 0.0)

    tx1, rx1 = LMP.modeterms(modeequation, tx, sampler)
    tx2, rx2 = LMP.modeterms(modeequation, tx, groundsampler)  # specialized

    @test tx1 ≈ tx2
    @test rx1 ≈ rx2

    @test length(tx1) == 1
    @test length(rx1) == LMP.NUMFIELDCOMPONENTS

    # frequency mismatch with modeequation
    txwrong = Transmitter(15e3)
    @test_throws ArgumentError LMP.modeterms(modeequation, txwrong, sampler)
end

function test_fieldsum(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    modes = TEST_MODES[scenario]

    X = LMP.distance(rx, tx)

    E1 = LMP.fieldsum(modes, waveguide, tx, rx)

    singlerx = GroundSampler(1000e3, rx.fieldcomponent)
    E2 = LMP.fieldsum(modes, waveguide, tx, singlerx)  # specialized

    distidx = findfirst(x->x==1000e3, X)
    @test E2 == E1[:,distidx]

    singlerx3 = GroundSampler(1000e3, Fields.E)
    E3 = LMP.fieldsum(modes, waveguide, tx, singlerx3)
    @test E3 == E2
end


@testset "modesum.jl" begin
    @info "Testing modesum"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        multiplespecies_scenario)
        
        test_modeterms(scn)
        test_fieldsum(scn)
    end
end
