using Test
using StaticArrays
using PolynomialRoots

using LongwaveModeSolver
const LWMS = LongwaveModeSolver


@testset "Integration Through Ionosphere" begin
    # Values taken from the end of a random homogeneous exponential LWPC run
    θ = complex(65.2520447, -1.5052794)
    ω = 150796.4531250  # rad/s

    Bfield = 0.5172267e-4  # (T) LWPC is actually using Gauss, despite their claim of W/m²=T
    dcl = -0.2664399
    dcm = -0.2850476
    dcn = -0.9207376

    referenceheight = 50.0
    topheight = 90.3906250
    bottomheight = 44.2968750
    height = 44.2968826

    @testset "Initial Condition (sharplyboundedreflectionmatrix)" begin
        M = @SMatrix [complex(-5.3366618, -8.5360718) complex(-5.8030024, -12.9980526) complex(-18.7323227, -27.8780575);
                      complex(-5.7961845, -5.0272179) complex(-6.1202874, -9.7538357) complex(-20.0426693, -32.2982941);
                      complex(-18.7344322, -30.3457165) complex(-20.0406971, -29.9917221) complex(-64.6526794, -100.7138519)]
        q = [complex(2.2052889, -0.0251389), complex(-0.0897404,-1.9890463)]
        r = @SMatrix [complex(2.7055206, -0.9693921) complex(0.7745517, 0.1726083);
                      complex(0.8630139, 0.2632782) complex(0.5952330, 0.6489610)]
        b = [1394.5465767833666 + 2201.647539088521im, 5.452661999120116 + 10.53265706490195im,
             119.74610270719461 + 133.6782356860242im, -34.67791183941958 - 52.48252736427487im,
             -63.6526794 - 100.7138519im]

        @test cosd(θ) ≈ complex(0.4187719, 0.0238619) atol=1e-6
        @test sind(θ) ≈ complex(0.9084715, -0.0109995) atol=1e-6

        qtest, Xtest = LWMS.sharplyboundedreflectionmatrix(θ, M)
        @test qtest[1:2] ≈ q atol=1e-6
        @test Xtest ≈ r atol=1e-6
    end

    @testset "M Matrix" begin
        electrons = Constituent(-fundamentalcharge, mₑ,
                                h -> 1.8548467e6, h -> 2.3626550e8)

        # 0s are not explicitly computed in LWPC
        Mlwpc = [complex(-0.0017909,-1.6545033E-004) 0 complex(1.9112933E-009,1.7586462E-006);
                 0 complex(-0.0017909,-1.6545285E-004) complex(-1.9975925E-009,-1.7646652E-006);
                 complex(-2.1218802E-009,-1.8791933E-006) complex(1.7722986E-009,1.6356993E-006) complex(-0.0017909,-1.6564084E-004)]

        Mtest, Ttest = LWMS.tmatrix(θ, ω, referenceheight, height, electrons, Bfield, dcl, dcm, dcn)
        @test Mtest[1,1] ≈ Mlwpc[1,1] atol=1e-7
        @test Mtest[2,2] ≈ Mlwpc[2,2] atol=1e-7
        @test Mtest[3,3] ≈ Mlwpc[3,3] atol=1e-7
        @test Mtest[1,3] ≈ Mlwpc[1,3] atol=1e-7
        @test Mtest[2,3] ≈ Mlwpc[2,3] atol=1e-7
        @test Mtest[3,1] ≈ Mlwpc[3,1] atol=1e-7

        @testset "S Matrix" begin
            A = [complex(3.9928365,-0.0006618) 0;
                 0 complex(4.0000000,0.0000000)]
            B = [complex(-0.8360517,-0.0475031) complex(-1.5016205E-008,-1.1766861E-005);
                 0 complex(-0.8375437,-0.0477239)]
            C = [complex(-0.8360517,-0.0475033) 0;
                 complex(-1.4886049E-008,-1.1692356E-005) complex(-0.8375437,-0.0477239)]
            D = [complex(0.0011740,3.7995025E-005) complex(-1.1764401E-007,3.9499610E-006);
                 complex(-1.5760942E-007,8.4526630E-007) complex(0.0017909,1.6545285E-004)]
            Atest, Btest, Ctest, Dtest = LWMS.smatrix(θ, Ttest)

            @test Atest ≈ A atol=1e-6
            @test Btest ≈ B atol=1e-6
            @test Ctest ≈ C atol=1e-6
            @test Dtest ≈ D atol=1e-6
        end
    end

    @testset "dXdh" begin
        X = [complex(2.7849803,-0.4330450) complex(0.1893668,0.3993837);
             complex(0.2346003,0.5070807) complex(2.4301479,0.2463614)]
        A = [complex(3.9928365,-0.0006618) 0;
             0 complex(4.0000000,0.0000000)]
        B = [complex(-0.8360517,-0.0475031) complex(-1.5016205E-008,-1.1766861E-005);
             0 complex(-0.8375437,-0.0477239)]
        C = [complex(-0.8360517,-0.0475033) 0;
             complex(-1.4886049E-008,-1.1692356E-005) complex(-0.8375437,-0.0477239)]
        D = [complex(0.0011740,3.7995025E-005) complex(-1.1764401E-007,3.9499610E-006);
             complex(-1.5760942E-007,8.4526630E-007) complex(0.0017909,1.6545285E-004)]
        k = ω/speedoflight*1e3

        lwpcderiv = [complex(0.1148383,0.1751488) complex(-0.1718566,0.0698231);
                     complex(-0.2180655,0.0862091) complex(-0.1612875,0.0093370)]

        testderiv = LWMS.dXdh(X, A, B, C, D, k)

        @test testderiv ≈ lwpcderiv atol=1e-6
    end

    @testset "Integration" begin
        electrons = Constituent(-fundamentalcharge, mₑ,
                                h -> waitprofile(h, 75, 0.3), collisionprofile)
        inputs = Inputs()
        inputs.topheight = topheight
        inputs.bottomheight = bottomheight

        Rbottomlwpc = [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
                       complex(0.2345901,0.5070829)  complex(2.4297192,0.2463784)]

        k = ω/speedoflight*1e3  # k must be based in km because we integrate over height in km

        sol = LWMS.integratethroughionosphere(θ, ω, k, topheight, bottomheight, referenceheight,
                                                    electrons, Bfield, dcl, dcm, dcn)
        @test sol[end] ≈ Rbottomlwpc atol=1e-1
    end
end

@testset "Free space integration" begin
    θ = complex(65.2520447,-1.5052794)
    ω = 150796.4531250  # rad/s
    k = ω/speedoflight*1e3

    referenceheight = 50.0
    fromheight = 44.2968750
    toheight = 67.3437500
    X₀ = [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
          complex(0.2345901,0.5070829) complex(2.4297192,0.2463784)]
    Xt = [complex(2.0947318,-0.0449929) complex(-0.0275066,-0.2524883);
          complex(-0.0319887,-0.3217891) complex(2.4113691,-0.3457983)]

    Xttest = LWMS.integratethroughfreespace!(X₀, θ, k, fromheight, toheight, referenceheight)

    @test Xttest ≈ Xt atol=1e-5
end

@testset "Ground reflection" begin
    σ = 0.003
    ϵᵣ = 15.
    ω = 150796.453125
    k = 0.5030022
    reflectionheight = 67.3437500
    referenceheight = 50.
    θ = complex(65.2520447,-1.5052794)

    num11 = complex(0.0423198,-0.2102124)
    num22 = complex(-0.0373752,0.8087913)
    den11 = complex(0.0210126,-0.1111768)
    den22 = complex(-0.0435943,0.2842564)

    N11, N22, D11, D22 = LWMS.groundreflection(θ, ω, k, σ, ϵᵣ, reflectionheight, referenceheight)

    @test N11 ≈ num11 atol=1e-5
    @test N22 ≈ num22 atol=1e-5
    @test D11 ≈ den11 atol=1e-5
    @test D22 ≈ den22 atol=1e-5
end

@testset "modifiedmodalfunction" begin
    σ = 0.003
    ϵᵣ = 15.
    ω = 150796.453125
    k = 0.5030022
    topheight = 90.3906250
    bottomheight = 44.2968750
    reflectionheight = 67.3437500
    referenceheight = 50.
    θ = complex(65.5882263,-0.9622505)

    B = 0.5172267e-4
    dcl = -0.2664398716543888
    dcm = -0.28504757220148297
    dcn = -0.9207375719361262

    electrons = Constituent(-fundamentalcharge, mₑ,
                            h -> waitprofile(h, 75, 0.3), collisionprofile)

    f = LWMS.modifiedmodalfunction(θ, ω, k, σ, ϵᵣ,
                              bottomheight, topheight, reflectionheight, referenceheight,
                              electrons,
                              B, dcl, dcm, dcn)

end
