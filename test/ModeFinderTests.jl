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

@testset "Free space integration" begin
    θ = complex(65.2520447,-1.5052794)
    ea = LWMS.EigenAngle(θ)
    ω = 150796.4531250  # rad/s
    freq = ω/2π
    source = LWMS.Source(freq)
    k = 1000source.k

    referenceheight = 50.0
    fromheight = 44.2968750
    toheight = 67.3437500
    X₀ = @MMatrix [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
                   complex(0.2345901,0.5070829) complex(2.4297192,0.2463784)]
    dXdC₀ = @MMatrix [complex(-15.3858480,-9.2181177) complex(11.8751659,-7.2130427);
                      complex(16.8020439,-9.8205147) complex(7.4359517,0.7653041)]
    Xt = [complex(2.0947318,-0.0449929) complex(-0.0275066,-0.2524883);
          complex(-0.0319887,-0.3217891) complex(2.4113691,-0.3457983)]
    dXt = [complex(-4.1203794,1.4973412) complex(-1.9841124,1.0458319);
           complex(-3.6317303,1.4457128) complex(-7.6892381,-0.8168259)]

    Xttest = LWMS.integratethroughfreespace!(X₀, ea, k, fromheight, toheight, referenceheight)

    @test Xttest ≈ Xt atol=1e-5

    X₀ = @MMatrix [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
                   complex(0.2345901,0.5070829) complex(2.4297192,0.2463784)]
    Xttest, dXttest = LWMS.integratethroughfreespace_XdXdC!(X₀, dXdC₀, ea, k, fromheight, toheight, referenceheight)

    @test Xttest ≈ Xt atol=1e-5
    @test dXttest ≈ dXt atol=1e-4
end

@testset "Ground reflection" begin
    σ = 0.003
    ϵᵣ = 15.
    ω = 150796.453125
    # k = 0.5030022
    freq = ω/2π
    source = LWMS.Source(freq)
    k = 1000source.k
    reflectionheight = 67.3437500
    referenceheight = 50.
    θ = complex(65.2520447,-1.5052794)
    ea = LWMS.EigenAngle(θ)

    num11 = complex(0.0423198,-0.2102124)
    num22 = complex(-0.0373752,0.8087913)
    den11 = complex(0.0210126,-0.1111768)
    den22 = complex(-0.0435943,0.2842564)

    groundcoeffs = LWMS._groundcoeffs(ea, source, σ, ϵᵣ, reflectionheight, referenceheight)
    N11, D11, N22, D22 = LWMS.ndground(groundcoeffs)

    @test N11 ≈ num11 atol=1e-5
    @test N22 ≈ num22 atol=1e-5
    @test D11 ≈ den11 atol=1e-5
    @test D22 ≈ den22 atol=1e-5

    C, C² = ea.cosθ, ea.cos²θ
    W, ng² = groundcoeffs.W, groundcoeffs.ng²
    n11, n22, d11, d22 = LWMS.fresnelnd(C, W, ng²)
    dn11dC, dd11dC, dn22dC, dd22dC = LWMS.dfresnelnddC(C, C², W, ng²)

    @test dn11dC ≈ complex(0.)
    @test dn22dC ≈ complex(2.9137998E-006,-2.6495843E-006) atol=1e-6
    @test dd11dC ≈ complex(0.5000014,-1.3253123E-006) atol=1e-6
    @test dd22dC ≈ complex(0.0074829,0.0074347) atol=1e-6

    dn11dC, dd11dC, dn22dC, dd22dC = LWMS.dndgrounddC(ea, groundcoeffs)

    @test dn11dC ≈ complex(-11.0291405,-1.0116535) atol=1e-3
    @test dn22dC ≈ complex(18.8285923,1.2034768) atol=1e-3
    @test dd11dC ≈ complex(-3.7694170,-1.0067097) atol=1e-3
    @test dd22dC ≈ complex(9.8127975,1.4418024) atol=1e-3
end
