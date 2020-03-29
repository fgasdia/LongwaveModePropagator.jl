using Test

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

#
# Profiles
#
h′ = 75
β = 0.32
N(z, h′, β) = 1.43e13*exp(-0.15*h′)*exp((β - 0.15)*(z/1000 - h′))

@test LWMS.waitprofile(80e3, h′, β) ≈ N(80e3, h′, β)


#
# BField
#
B = 50000e-9
dip = deg2rad(-78)
az = deg2rad(240)
bfield = LWMS.BField(B, dip, az)

@test LWMS.dip(bfield) ≈ dip
@test LWMS.azimuth(bfield) ≈ az

dip = deg2rad(62.4)
az = deg2rad(18.3)
bfield = LWMS.BField(B, dip, az)
@test LWMS.dip(bfield) ≈ dip
@test LWMS.azimuth(bfield) ≈ az

# TODO: Test if direction cosines are correct
@test_skip bfield.dcl == 2
