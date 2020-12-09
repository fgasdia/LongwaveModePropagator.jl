function tmatrix_deriv(scenario)
    @unpack ea, tx, bfield, species = scenario()

    M = LMP.susceptibility(LMPParams().topheight, tx.frequency, bfield, species)

    for i = 1:4
        for j = 1:4
            Tfcn(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M)[i,j])
            dTref = FiniteDiff.finite_difference_derivative(Tfcn, θs, Val{:central})
            dT(θ) = (ea = EigenAngle(θ); T = LMP.tmatrix(ea, M, LMP.Dθ())[i,j])

            @test err_func(dT.(θs), dTref) < 1e-6
        end
    end
end

@testset "TMatrix.jl" begin
    @info "Testing TMatrix"

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario)
        tmatrix_deriv(scn)
    end
end
