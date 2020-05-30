function ne_waitprofile()
    hps = 60.0:90.0
    betas = 0.15:0.05:2
    zs = 10e3:2e3:90e3

    N(z, h′, β) = 1.43e13*exp(-0.15*h′)*exp((β - 0.15)*(z/1000 - h′))

    for h in hps, β in betas, z in zs
        waitprofile(z, h, β) ≈ N(z, h, β) || return false
    end

    # Check if zero below cutoff height
    waitprofile(30e3, 60, 0.35, 50e3) == 0 || return false

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

@testset "Geophysics" begin
    @info "Testing geophysics functions..."

    @testset "Profiles" begin
        @info "  Testing profiles..."

        @test ne_waitprofile()
    end

    @testset "Magnetic field" begin
        @info "  Testing magnetic fields..."

        @test magnetic_dipaz()
        @test magnetic_directioncosines()
    end
end
