using Test
using StaticArrays
using PolynomialRoots

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

    referenceheight = 50.0
    topheight = 90.3906250
    bottomheight = 44.2968750
    height = 44.2968826

    X = @MMatrix zeros(ComplexF64, 2, 2)
    dXdC = similar(X)

    @testset "Initial Condition (sharplybounded_X)" begin
        M = @SMatrix [complex(-5.3366618, -8.5360718) complex(-5.8030024, -12.9980526) complex(-18.7323227, -27.8780575);
                      complex(-5.7961845, -5.0272179) complex(-6.1202874, -9.7538357) complex(-20.0426693, -32.2982941);
                      complex(-18.7344322, -30.3457165) complex(-20.0406971, -29.9917221) complex(-64.6526794, -100.7138519)]
        q = [complex(2.2052889, -0.0251389), complex(-0.0897404,-1.9890463)]
        r = @SMatrix [complex(2.7055206, -0.9693921) complex(0.7745517, 0.1726083);
                      complex(0.8630139, 0.2632782) complex(0.5952330, 0.6489610)]
        b = [1394.5465767833666 + 2201.647539088521im, 5.452661999120116 + 10.53265706490195im,
             119.74610270719461 + 133.6782356860242im, -34.67791183941958 - 52.48252736427487im,
             -63.6526794 - 100.7138519im]

        qtest, bcoeffs = LWMS.bookerquartic(ea, M)
        D = LWMS._sharplybounded(ea, M, q)
        Xtest = LWMS.sharplyboundedX!(X, D, ea, q)
        @test qtest[1:2] ≈ q atol=1e-6
        @test Xtest ≈ r atol=1e-6

        dr11 = complex(-3.4461095,2.2281921)
        dr22 = complex(-0.1982546,-0.5395676)
        dr12 = complex(-1.1654160,-0.2949337)
        dr21 = complex(-1.4740430,-0.3747938)

        dXtest = LWMS.sharplyboundeddXdC!(dXdC, Xtest, D, ea, M, q, bcoeffs)

        # LWPC likely has a bug in the derivative
        @test dXtest[1,1] ≈ dr11 atol=1
        @test dXtest[2,1] ≈ dr21 atol=1e-1
        @test dXtest[1,2] ≈ dr12 atol=1e-1
        @test dXtest[2,2] ≈ dr22 atol=1e-1
    end

    @testset "M Matrix" begin
        electrons = Constituent(-fundamentalcharge, mₑ,
                                h -> 1.8548467e6, h -> 2.3626550e8)

        # 0s are not explicitly computed in LWPC
        Mlwpc = [complex(-0.0017909,-1.6545033E-004) 0 complex(1.9112933E-009,1.7586462E-006);
                 0 complex(-0.0017909,-1.6545285E-004) complex(-1.9975925E-009,-1.7646652E-006);
                 complex(-2.1218802E-009,-1.8791933E-006) complex(1.7722986E-009,1.6356993E-006) complex(-0.0017909,-1.6564084E-004)]

        Mtest = @MMatrix zeros(ComplexF64, 3, 3)

        Mtest = LWMS.mmatrix!(Mtest, ω, referenceheight, height, electrons, bfield)
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
            Atest, Btest, Ctest, Dtest = LWMS.smatrix(ea, Mtest)

            @test Atest ≈ A atol=1e-6
            @test Btest ≈ B atol=1e-6
            @test Ctest ≈ C atol=1e-6
            @test Dtest ≈ D atol=1e-6

            dB = [complex(-1.9964184,3.3252075E-004) 0;
                  0 complex(-2.0000000,0.0000000)]
            dC = [complex(-1.9964184,3.3263181E-004) 0;
                  0 complex(-2.0000000,0.0000000)]
            dD = [complex(-0.0029866,-4.4628372E-004) complex(5.8776678E-008,5.1300080E-006);
                  complex(-4.7828390E-008,6.6590064E-006) 0]

            dBtest, dCtest, dDtest = LWMS.dsmatrixdC(ea, Mtest)
            @test dBtest ≈ dB atol=1e-6
            @test dCtest ≈ dC atol=1e-6
            @test dDtest ≈ dD atol=1e-6
        end
    end

    @testset "dXdh" begin
        X = @SMatrix [complex(2.7849803,-0.4330450) complex(0.1893668,0.3993837);
             complex(0.2346003,0.5070807) complex(2.4301479,0.2463614)]
        A = @SMatrix [complex(3.9928365,-0.0006618) 0;
             0 complex(4.0000000,0.0000000)]
        B = @SMatrix [complex(-0.8360517,-0.0475031) complex(-1.5016205E-008,-1.1766861E-005);
             0 complex(-0.8375437,-0.0477239)]
        C = @SMatrix [complex(-0.8360517,-0.0475033) 0;
             complex(-1.4886049E-008,-1.1692356E-005) complex(-0.8375437,-0.0477239)]
        D = @SMatrix [complex(0.0011740,3.7995025E-005) complex(-1.1764401E-007,3.9499610E-006);
             complex(-1.5760942E-007,8.4526630E-007) complex(0.0017909,1.6545285E-004)]

        k = 1000source.k

        lwpcderiv = [complex(0.1148383,0.1751488) complex(-0.1718566,0.0698231);
                     complex(-0.2180655,0.0862091) complex(-0.1612875,0.0093370)]

        testderiv = LWMS.dXdh(X, A, B, C, D, k)

        @test testderiv ≈ lwpcderiv atol=1e-6

        dB = @SMatrix [complex(-1.9964184,3.3252075E-004) 0;
                        0   complex(-2.0000000,0.0000000)]
        dC = @SMatrix [complex(-1.9964184,3.3263181E-004) 0;
                        0 complex(-2.0000000,0.0000000)]
        dD = @SMatrix [complex(-0.0029866,-4.4628372E-004) complex(5.8776678E-008,5.1300080E-006);
                        complex(-4.7828390E-008,6.6590064E-006) 0]

        dX = @SMatrix [complex(-15.3878374,-9.2189331) complex(11.8746481,-7.2133737);
                       complex(16.8013687,-9.8209696) complex(7.4344139,0.7641335)]

        lwpcderivc = [complex(4.6734123,-3.4267042) complex(2.3361709,5.3375201);
                      complex(3.2026398,7.5093975) complex(-0.7406456,5.5385509)]

        testderivc = LWMS.dXdCdh(X, dX, B, dB, C, dC, D, dD, k)

        @test testderivc ≈ lwpcderivc atol=1e-4
    end

    @testset "Integration" begin
        electrons = Constituent(-fundamentalcharge, mₑ,
                                h -> waitprofile(h, 75, 0.3), collisionprofile)
        inputs = Inputs()
        inputs.topheight = topheight
        inputs.bottomheight = bottomheight

        Rbottomlwpc = [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
                       complex(0.2345901,0.5070829)  complex(2.4297192,0.2463784)]

        sol = LWMS.integratethroughionosphere(ea, source, topheight, bottomheight, referenceheight,
                                              electrons, bfield)
        @test sol[end] ≈ Rbottomlwpc atol=1e-1
    end

    @testset "Integration with dC" begin
        electrons = Constituent(-fundamentalcharge, mₑ,
                                h -> waitprofile(h, 75, 0.3), collisionprofile)
        inputs = Inputs()
        inputs.topheight = topheight
        inputs.bottomheight = bottomheight

        Rbottomlwpc = [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
                       complex(0.2345901,0.5070829)  complex(2.4297192,0.2463784)]

        sol = LWMS.integratethroughionosphere_dC(ea, source, topheight, bottomheight, referenceheight,
                                                 electrons, bfield)
        Xsol = @view sol[end][1:2,:]
        dXsol = @view sol[end][3:4,:]

        dXlwpc = [complex(-15.3858480,-9.2181177) complex(11.8751659,-7.2130427);
                  complex(16.8020439,-9.8205147) complex(7.4359517,0.7653041)]

        @test Xsol ≈ Rbottomlwpc atol=1e-1
        @test dXsol ≈ dXlwpc atol=1  # LWPC has a bug in calculation of init dX?
    end
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
    Xt = [complex(2.0947318,-0.0449929) complex(-0.0275066,-0.2524883);
          complex(-0.0319887,-0.3217891) complex(2.4113691,-0.3457983)]

    Xttest = LWMS.integratethroughfreespace!(X₀, ea, k, fromheight, toheight, referenceheight)

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

    groundcoeffs = LWMS._groundcoeffs(ea, source, σ, ϵᵣ, reflectionheight, referenceheight)
    N11, N22, D11, D22 = LWMS.ndground(groundcoeffs)

    @test N11 ≈ num11 atol=1e-5
    @test N22 ≈ num22 atol=1e-5
    @test D11 ≈ den11 atol=1e-5
    @test D22 ≈ den22 atol=1e-5

    C = ea.cosθ
    W, ng² = groundcoeffs.W, groundcoeffs.ng²
    n11, n22, d11, d22 = LWMS.fresnelnd(C, W, ng²)
    dn11dC, dn22dC, dd11dC, dd22dC = LWMS.dfresnelnddC(C, W, ng²)

    @test dn11dC ≈ complex(0.)
    @test dn22dC ≈ complex(2.9137998E-006,-2.6495843E-006) atol=1e-6
    @test dd11dC ≈ complex(0.5000014,-1.3253123E-006) atol=1e-6
    @test dd22dC ≈ complex(0.0074829,0.0074347) atol=1e-6

    dn11dC, dn22dC, dd11dC, dd22dC = LWMS.dndgrounddC(ea, groundcoeffs)

    @test dn11dC ≈ complex(-11.0291405,-1.0116535) atol=1e-3
    @test dn22dC ≈ complex(18.8285923,1.2034768) atol=1e-3
    @test dd11dC ≈ complex(-3.7694170,-1.0067097) atol=1e-3
    @test dd22dC ≈ complex(9.8127975,1.4418024) atol=1e-3
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
                              bfield)

end
