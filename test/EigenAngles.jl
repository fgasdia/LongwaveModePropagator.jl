@testset "EigenAngles.jl" begin
    @info "Testing EigenAngles"

    @test EigenAngle(deg2rad(80-0.5im)) > EigenAngle(deg2rad(75-0.3im))

    @test_logs (:warn, "θ > 2π. Make sure θ is in radians.") EigenAngle(complex(85.0, 0.31))
end
