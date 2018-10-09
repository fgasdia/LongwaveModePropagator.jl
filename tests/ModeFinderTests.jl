"""
ModeFinder tests
"""
module ModeFinderTests

# thisdir = dirname(@__FILE__())
# any(path -> path==thisdir, LOAD_PATH) || push!(LOAD_PATH, thisdir)
# srcdir = abspath(thisdir, "..\\src")
# any(path -> path==srcdir, LOAD_PATH) || push!(LOAD_PATH, srcdir)

include("../src/ModeFinder.jl")
include("../src/LWMS.jl")

using Test

using StaticArrays
using PolynomialRoots

using LWMS
using ModeFinder
using IonosphereProfile


"""
Solution for the roots of a fourth-order polynomial (quartic equation) taken from Burnside
and Panton (1904), the Theory of Equations. A summary of the equations and process is given
in Radio Science Vol. 3, Aug 1968, pg. 792-795.
"""
function solvequartic(B4, B3, B2, B1, B0)
    # More suitable form for solution:
    # q⁴ + 4b₃q³ + 6b₂q² + 4b₁q + b₀ = 0
    b3 = B3/4B4
    b2 = B2/6B4
    b1 = B1/4B4
    b0 = B0/B4

    H = b2 - b3^2
    I = b0 - 4*b3*b1 + 3b2^2
    G = b1 - 3*b3*b2 + 2b3^3
    h = -I/12
    g = -G^2/4 - H*(H^2 + 3h)
    p = (-g + sqrt(g^2 + 4h^3))/2
    magpos = abs(real(p)) + abs(imag(p))
    ppos = p
    p = (-g - sqrt(g^2 + 4h^3))/2
    magneg = abs(real(p)) + abs(imag(p))
    magpos > magneg && (p = ppos)

    ω2 = complex(-0.5, sqrt(3)/2)
    ω3 = complex(-0.5, -sqrt(3)/2)
    s1 = exp(log(p)/3)  # cube root of p
    s2 = ω2*s1
    s3 = ω3*s1

    r1 = sqrt(s1 - h/s1 - H)
    r2 = sqrt(s2 - h/s2 - H)
    r3 = sqrt(s3 - h/s3 - H)

    if abs(G) >= 1e-20
        real(-2r1*r2*r3/G) ≈ -1 && (r3 = -r3)
    end

    q = @MVector [r1+r2+r3, r1-r2-r3, -r1+r2-r3, -r1-r2+r3]
    q .-= b3

    # First order Newton-Raphson iterative improvement
    for j = 1:4
        dlqmin = 9.9e9
        for jj = 1:4
            if jj != j
                dlq = abs(q[j] - q[jj])
                dlq < dlqmin && (dqlmin = dlq)
            end
        end
        dlqmax = dlqmin/3

        lastiter = false
        ncount = 1
        delq = 0.0
        while !lastiter & (ncount <= 10)
            f = (((B4*q[j] + B3)*q[j] + B2)*q[j] + B1)*q[j] + B0
            dfdq = ((4B4*q[j] + 3B3)*q[j] + 2B2)*q[j] + B1
            delq = -f/dfdq
            abs(delq) > dlqmax && (delq *= delqmax/abs(delq))
            q[j] += delq
            abs(delq/q[j]) < 1e-4 && (lastiter = true)
            ncount += 1
        end
        if abs(delq/q[j]) > 1e-2
            error("q fails to converge")
        end
    end
    return q
end

@testset "Integration Through Ionosphere" begin
    # Values taken from the end of a random homogeneous exponential LWPC run
    θ = complex(65.2520447, -1.5052794)
    ω = 150796.4531250  # rad/s

    Bfield = 0.5172267e-4  # (T) LWPC is actually using Gauss, despite their claim of W/m²=T
    dcl = -0.2664399
    dcm = -0.2850476
    dcn = -0.9207376

    ht = 44.2968826
    topheight = 90.3906250
    bottomheight = 44.2968750
    referenceheight = 50.0


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

        qtest, Xtest = ModeFinder.sharplyboundedreflectionmatrix(θ, M)
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

        Mtest, Ttest = ModeFinder.tmatrix(θ, ht, referenceheight, ω, electrons, Bfield, dcl, dcm, dcn)
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
            Atest, Btest, Ctest, Dtest = ModeFinder.smatrix(θ, Ttest)

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

        testderiv = ModeFinder.dXdh(X, A, B, C, D, k)

        @test testderiv ≈ lwpcderiv atol=1e-6
    end

    @testset "Integration" begin
        electrons = Constituent(-fundamentalcharge, mₑ,
                                h -> waitprofile(h, 75, 0.3), collisionprofile)
        inputs = Inputs()
        inputs.topheight = topheight
        inputs.bottomheight = bottomheight
        drcs = ModeFinder.DirectionCosines(dcl, dcm, dcn, 0.0)

        Rbottomlwpc = [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
                       complex(0.2345901,0.5070829)  complex(2.4297192,0.2463784)]

        mode = ModeFinder.Mode()
        mode.θ = θ
        mode.ω = ω
        mode.wavenumber = ω/speedoflight*1e3  # k must be based in km because we integrate over height in km

        sol = ModeFinder.integratethroughionosphere(mode, inputs, referenceheight,
                                                    electrons, Bfield, drcs)
        @test sol[end] ≈ Rbottomlwpc atol=1e-1
    end
end

@testset "Free space integration" begin
    @testset "Modified Hankel" begin
        z = complex(24.6770630,2.7361517)
        h1 = complex(0.3822204,-0.0108716)
        h2 = complex(0.3823441,-0.0102396)
        h1p = complex(-0.0548274,1.9051478)
        h2p = complex(0.0503774,-1.9045298)

        # h1test = ModeFinder.modhankel1(z)
        # h2test = ModeFinder.modhankel2(z)
        # h1ptest = ModeFinder.modhankel1prime(z)
        # h2ptest = ModeFinder.modhankel2prime(z)

        #  Note: this function is not currently in use
        mh1test, mh2test, mh1ptest, mh2ptest = ModeFinder.modhankel(z, false)

        @test mh1test ≈ h1 atol=1e-6
        @test mh2test ≈ h2 atol=1e-6
        @test mh1ptest ≈ h1p atol=1e-6
        @test mh2ptest ≈ h2p atol=1e-6
    end

    mode = ModeFinder.Mode()
    mode.θ = complex(65.2520447,-1.5052794)
    mode.ω = 150796.4531250  # rad/s
    mode.wavenumber = mode.ω/speedoflight*1e3

    h = 50.0
    z₀ = 44.2968750
    zz = 67.3437500
    X₀ = [complex(2.7846940,-0.4330095) complex(0.1893602,0.3993874);
          complex(0.2345901,0.5070829) complex(2.4297192,0.2463784)]
    Xt = [complex(2.0947318,-0.0449929) complex(-0.0275066,-0.2524883);
          complex(-0.0319887,-0.3217891) complex(2.4113691,-0.3457983)]

    Xttest = ModeFinder.integratethroughfreespace!(X₀, mode, z₀, zz, h)

    @test Xttest ≈ Xt atol=1e-2
end


@testset "Ordinary/Extraordinary Waves" begin
    @testset "Eigenmatrix" begin
        dcl = -0.2664399
        dcm = -0.2850476
        dcn = -0.9207376
        G = complex(-0.1100280,2.0067947E-004)
        h = 50.0
        ec_lwpc = 4.0721876E-004
        heightinterp = ec_lwpc*ModeFinder.earthradius/2+h

        lenset = false

        mode = ModeFinder.Mode()
        mode.ω = 150796.4531250  # rad/s
        mode.θ = complex(39.9019470,-1.7546921)

        drcs = ModeFinder.DirectionCosines(dcl, dcm, dcn, G)

        Eud = @MArray zeros(ComplexF64, 2,2,2)

        ζₚ = @MVector [complex(0.0389593,-1.7352930), complex(0.0773762,-1.0256317)]
        Γₚ = @MVector [complex(2.4567194,0.0617389), complex(1.4542819,0.1539752)]
        xseq = ModeFinder.XSequence(ζₚ, Γₚ)

        e_lwpc = Array{ComplexF64}(undef,2,2,2)
        e_lwpc[1,1,1] = complex(-0.6970944,0.0006473)
        e_lwpc[2,1,1] = complex(0.0014815,-0.7075908)
        e_lwpc[1,2,1] = complex(0.0014815,-0.7075908)
        e_lwpc[2,2,1] = complex(-0.7162860,0.0023426)
        e_lwpc[1,1,2] = complex(0.7391026,-0.0233032)
        e_lwpc[2,1,2] = complex(0.0199049,-0.7092856)
        e_lwpc[1,2,2] = complex(0.0199049,-0.7092856)
        e_lwpc[2,2,2] = complex(0.6733917,-0.0169724)

        ModeFinder.eigenmatrix!(Eud, xseq, mode, lenset, heightinterp, h, drcs)

        @test Eud ≈ e_lwpc atol=1e-6
    end

    @testset "ROE/R matrix conversions" begin
        θ = complex(42.9411621,-3.6225901)
        X = [complex(1.7967621,-0.0265428) complex(0.0993883,0.3174350);
             complex(0.5336750,0.3858191) complex(1.2834927,0.4988435)]

        Eud = zeros(ComplexF64, 2,2,2)
        Eud[1,1,1] = complex(-0.6950788,0.0005265)
        Eud[2,1,1] = complex(0.0023001,-0.7082505)
        Eud[1,2,1] = complex(0.0023001,-0.7082505)
        Eud[2,2,1] = complex(-0.7170246,0.0041443)
        Eud[1,1,2] = complex(0.7401037,-0.0300393)
        Eud[2,1,2] = complex(0.0229389,-0.7119964)
        Eud[1,2,2] = complex(0.0229389,-0.7119964)
        Eud[2,2,2] = complex(0.6676079,-0.0170387)

        R_lwpc = [complex(0.3191003,0.0579732) complex(0.0592212,0.2371279);
                  complex(0.3748306,0.3060071) complex(-0.0800366,0.4212306)]
        Roe_lwpc = [complex(0.0721194,-0.0232919) complex(0.3896856,-0.0565601);
                    complex(0.0823480,-0.1521263) complex(0.4669252,-0.4200467)]

        Rtest, Roetest = ModeFinder.roematrix(θ, X, Eud)
        @test Roetest ≈ Roe_lwpc atol=1e-6
        @test Rtest ≈ R_lwpc atol=1e-6

        Rtest, Xtest = ModeFinder.rmatrix(θ, Roe_lwpc, Eud)
        @test Xtest ≈ X atol=1e-6
        @test Rtest ≈ R_lwpc atol=1e-6
    end
end

end
