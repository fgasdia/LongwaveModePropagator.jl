function test_tmatrix_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Tfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M)[i,j])
            dTref = FiniteDiff.finite_difference_derivative(Tfcn, θs, Val{:central})
            dT(θ) = (ea = EigenAngle(θ); T = LMP.dtmatrix(ea, M)[i,j])

            @test maxabsdiff(dT.(θs), dTref) < 1e-7
        end
    end
end

@testset "TMatrix.jl" begin
    @info "Testing TMatrix"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        multiplespecies_scenario)
        test_tmatrix_deriv(scn)
    end
end
