using Test
using StaticArrays

using LongwaveModeSolver
const LWMS = LongwaveModeSolver


@testset "Integration Through Ionosphere" begin
    # Values taken from the end of a random homogeneous exponential LWPC run
    θ = complex(65.2520447, -1.5052794)
    ω = 150796.4531250  # rad/s
    freq = ω/2π

    ea = LWMS.EigenAngle(θ)
    source = LWMS.Source(freq)

    Bmag = 0.5172267e-4  # (T) LWPC is actually using Gauss, despite their claim of W/m²=T
    dcl = -0.2664399
    dcm = -0.2850476
    dcn = -0.9207376
    bfield = LWMS.BField(Bmag, dcl, dcm, dcn)

    electrons = Constituent(-fundamentalcharge, mₑ,
                            h -> waitprofile(h, 75, 0.3), collisionprofile)

    referenceheight = 50.0
    topheight = 90.3906250
    bottomheight = 44.2968750
    height = 70.0

    sol = LWMS.integratethroughionosphere(ea, source, topheight, bottomheight, referenceheight,
                                          electrons, bfield)
end
