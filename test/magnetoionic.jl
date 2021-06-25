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

function test_spline(scenario)
    @unpack tx, bfield, species = scenario

    itp = LMP.susceptibilityspline(tx.frequency, bfield, species)

    zs = LMP.BOTTOMHEIGHT:1:LMPParams().topheight
    Ms = LMP.susceptibility.(zs, (tx.frequency,), (bfield,), (species,))
    Mitp = itp.(zs)

    function matrix(M)
        mat = Matrix{eltype(M[1])}(undef, length(M), 9)
        for i in eachindex(M)
            mat[i,:] = M[i]
        end
        return mat
    end

    @test matrix(Mitp) ≈ matrix(Ms) rtol=1e-5

    # plot(matrix(real(Ms)), zs/1000, legend=false)
    # plot!(matrix(real(Mitp)), zs/1000, linestyle=:dash)

    # Mdiff = matrix(real(Ms)) .- matrix(real(Mitp))
    # scaleddiff = Mdiff./matrix(real(Ms))
    # plot(Mdiff, zs/1000, legend=false, ylims=(0,110))

    # Make sure `params.susceptibilitysplinestep` is being used
    itp2 = LMP.susceptibilityspline(tx.frequency, bfield, species; params=LMPParams(susceptibilitysplinestep=10.0))
    itp3 = LMP.susceptibilityspline(tx.frequency, bfield, species; params=LMPParams(susceptibilitysplinestep=20.0))
    @test itp == itp2
    @test itp != itp3
end

@testset "magnetoionic.jl" begin
    @info "Testing magnetoionic"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
            multiplespecies_scenario)
        test_susceptibility(scn)
        test_spline(scn)
    end
end
