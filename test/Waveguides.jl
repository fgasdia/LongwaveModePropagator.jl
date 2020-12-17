function test_HomogeneousWaveguide(scenario)
    @unpack bfield, species, ground = scenario

    distance = 1000e3
    waveguide = HomogeneousWaveguide(bfield, species, ground, distance)
    @test isbits(waveguide)

    adjwaveguide = LMP.adjoint(waveguide)
    oB = waveguide.bfield
    aB = adjwaveguide.bfield
    @test aB.B == oB.B
    @test aB.dcl == -oB.dcl
    @test aB.dcm == oB.dcm
    @test aB.dcn == oB.dcn
    @test adjwaveguide.species == waveguide.species
    @test adjwaveguide.ground == waveguide.ground
    @test adjwaveguide.distance == waveguide.distance
end

@testset "Waveguides.jl" begin
    @info "Testing Waveguides"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        test_HomogeneousWaveguide(scn)
    end
end
