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
