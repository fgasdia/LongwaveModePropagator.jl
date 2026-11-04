function test_susceptibility(scenario)
    @unpack tx, bfield, species, ground = scenario

    M1 = LMP.susceptibility(70e3, tx.frequency, bfield, species)
    M2 = LMP.susceptibility(70e3, tx.frequency, bfield, species; params=LMPParams())
    M3 = LMP.susceptibility(70e3, tx.frequency, bfield, species; params=LMPParams(earthradius=6350e3))
    @test M1 == M2
    @test !(M2 ≈ M3)

    waveguide = HomogeneousWaveguide(bfield, species, ground)
    modeequation = PhysicalModeEquation(tx.frequency, waveguide)

    M4 = LMP.susceptibility(70e3, tx.frequency, waveguide)
    M5 = LMP.susceptibility(70e3, tx.frequency, waveguide; params=LMPParams())
    M6 = LMP.susceptibility(70e3, tx.frequency, waveguide; params=LMPParams(earthradius=6350e3))

    M7 = LMP.susceptibility(70e3, modeequation)
    M8 = LMP.susceptibility(70e3, modeequation; params=LMPParams())
    M9 = LMP.susceptibility(70e3, modeequation; params=LMPParams(earthradius=6350e3))

    @test M4 == M5 == M1
    @test M6 == M3

    @test M7 == M8 == M1
    @test M9 == M3

    @inferred LMP.susceptibility(70e3, tx.frequency, bfield, species)
end

@testset "magnetoionic.jl" begin
    @info "Testing magnetoionic"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
            multiplespecies_scenario)
        test_susceptibility(scn)
    end
end
