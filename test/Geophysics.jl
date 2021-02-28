function test_Species()
    electrons = Species(QE, ME,
                        z -> waitprofile(z, 75, 0.32), electroncollisionfrequency)
    @test isbits(electrons)
end

function test_BField()
    bfield = BField(50e-6, π/2, 3π/2)
    @test isbits(bfield)
    @test (@test_logs (:warn, "BField magnitude of exactly 0 is not supported."*
                       " Setting B = 1e-16.") BField(0, π/2, 3π/2).B) == 1e-16

    #==
    dip and az
    ==#
    B = 50000e-9
    true_dip = deg2rad(-78)
    true_az = deg2rad(240)
    bfield = BField(B, true_dip, true_az)

    @test LMP.dip(bfield) ≈ true_dip
    @test mod2pi(LMP.azimuth(bfield)) ≈ true_az

    true_dip = deg2rad(62.4)
    true_az = deg2rad(18.3)
    bfield = BField(B, true_dip, true_az)

    @test LMP.dip(bfield) ≈ true_dip
    @test mod2pi(LMP.azimuth(bfield)) ≈ true_az

    #==
    direction cosines
    ==#
    # dip angle 90° up from horizontal
    bfield = BField(50e-6, π/2, 3π/2)
    @test bfield.dcn == -1

    # dip angle 90° down from horizontal
    bfield = BField(50e-6, -π/2, 3π/2)
    @test bfield.dcn == 1
end

function test_Ground()
    @test isbits(Ground(15, 0.001))

    @test LMP.GROUND[3] == Ground(10, 1e-4)
end

function test_waitprofile()
    zs = 50e3:2e3:86e3
    hbs = [(65, 0.2), (65, 0.25), (68, 0.25), (70, 0.25), (72, 0.3), (72, 0.35), (75, 0.35),
           (78, 0.35), (78, 0.4), (82, 0.5), (82, 0.55), (85, 0.6), (85, 0.65), (88, 0.8)]

    N(z, h′, β) = 1.43e13*exp(-0.15*h′)*exp((β - 0.15)*(z/1000 - h′))

    for hb in hbs, z in zs
        @test waitprofile(z, hb[1], hb[2]) ≈ N(z, hb[1], hb[2])
    end

    # Check if zero below cutoff height
    @test waitprofile(30e3, 60, 0.35; cutoff_low=40e3) == 0
    @test waitprofile(30e3, 60, 0.35; cutoff_low=40e3) isa Float64
    @test waitprofile(30_000, 60, 0.35; cutoff_low=40e3) isa Float64

    # Check default threshold
    @test waitprofile(110e3, 70, 0.5) == 1e12
end

@testset "Geophysics.jl" begin
    @info "Testing Geophysics"

    test_Species()
    test_BField()
    test_Ground()
    test_waitprofile()
    @test electroncollisionfrequency(110e3) isa Float64
    @test ioncollisionfrequency(110e3) isa Float64
end
