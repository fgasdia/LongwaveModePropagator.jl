function bitstype()
    ea = EigenAngle(deg2rad(complex(85.0, -1.0)))
    return isbits(ea)
end

function warns_radians()
    EigenAngle(complex(85.0, 0.31))
end

function sorting()
    EigenAngle(deg2rad(80-0.5im)) > EigenAngle(deg2rad(75-0.3im))
end

@testset "EigenAngles.jl" begin
    @info "Testing EigenAngles"

    @test bitstype()
    @test sorting()

    @test_logs (:warn, "θ should be in radians") warns_radians()
end
