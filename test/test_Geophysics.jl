function ne_waitprofile()
    zs = 50e3:2e3:86e3
    hbs = [(65, 0.2), (65, 0.25), (68, 0.25), (70, 0.25), (72, 0.3), (72, 0.35), (75, 0.35),
           (78, 0.35), (78, 0.4), (82, 0.5), (82, 0.55), (85, 0.6), (85, 0.65), (88, 0.8)]

    N(z, h′, β) = 1.43e13*exp(-0.15*h′)*exp((β - 0.15)*(z/1000 - h′))

    for hb in hbs, z in zs
        waitprofile(z, hb[1], hb[2]) ≈ N(z, hb[1], hb[2]) || return false
    end

    return true
end

function magnetic_dipaz()
    B = 50000e-9
    true_dip = deg2rad(-78)
    true_az = deg2rad(240)
    bfield = BField(B, true_dip, true_az)

    dip(bfield) ≈ true_dip || return false
    mod2pi(azimuth(bfield)) ≈ true_az || return false

    true_dip = deg2rad(62.4)
    true_az = deg2rad(18.3)
    bfield = BField(B, true_dip, true_az)

    dip(bfield) ≈ true_dip || return false
    mod2pi(azimuth(bfield)) ≈ true_az || return false

    return true
end

function magnetic_directioncosines()
    # dip angle 90° up from horizontal
    bfield = BField(50e-6, π/2, 3π/2)

    bfield.dcn == -1 || return false

    # dip angle 90° down from horizontal
    bfield = BField(50e-6, -π/2, 3π/2)

    bfield.dcn == 1 || return false

    # TODO: more tests...

    return true
end

function speciesbits()
    electrons = Species(QE, ME,
                        z -> waitprofile(z, 75, 0.32),
                        z -> electroncollisionfrequency)
    return isbits(electrons)
end

@testset "Geophysics.jl" begin
    @info "Testing Geophysics"

    @testset "Profiles" begin
        @info "  Profiles..."

        @test ne_waitprofile()

        # Check if zero below cutoff height
        @test waitprofile(30e3, 60, 0.35, cutoff_low=40e3) == 0 || return false

        # Check default threshold
        @test waitprofile(110e3, 70, 0.5) == 1e12 || return false
    end

    @testset "Species" begin
        @info "  Species..."

        @test speciesbits()
    end

    @testset "Magnetic field" begin
        @info "  BField..."

        @test magnetic_dipaz()
        @test magnetic_directioncosines()
        @test isbits(BField(50e-6, π/2, 3π/2))
        @test BField(0, π/2, 3π/2).B == 1e-15
        @test_logs (:warn, "B field magnitude of exactly 0 is not supported. Setting B = 1e-15.") BField(0, π/2, 3π/2);
    end

    @testset "Ground" begin
        @info "  Ground..."

        @test isbits(Ground(15, 0.001))
    end
end
