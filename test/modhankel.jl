
"""
This function calculates the modified Hankel functions of order 1/3 and their derivatives.

The functions h₁(ζ) and h₂(ζ) satisfy the Stokes differential equation (_Airy_ function)
```math
\\frac{d²w}{dζ²} + ζw = 0
```

Note: This may actually be doing more than just calculating the hankel functions, because
when used with `espace()` below, which follows the math steps in MS 1976,
the values from this function, `modhankel()`, result in the wrong values even though this
function matches the LWPC routine mf_mdhnkl. However, when the above functions `modhankel1`,
etc are used in `integratethroughfreespace`, that function agrees with LWPC mf_fsinteg.

For some reason, the transformation from Airy Ai and Bi to Mod Hankel is breaking for large
`z`, so I need to use LWPC's asymptotic expansion technique in those conditions.

TODO: Also not sure what LWPC `ngboth` is..., although it is fully supported right now because it
is only used in the large `z` case. It might mean 'negative both' based on use in mf_fsinteg
"""
function modhankel(z, bothnegative::Bool)
    # TODO: Why is the second part necessary? The airy function should already contain an
    # asymptotic expansion

    # XXX: No idea where this comes from
    cap = [1.0416666666666666663e-01,  8.3550347222222222116e-02,
           1.2822657455632716019e-01,  2.9184902646414046315e-01,
           8.8162726744375764874e-01,  3.3214082818627675264e+00,
           1.4995762986862554546e+01,  7.8923013011586517530e+01,
           4.7445153886826431887e+02,  3.2074900908906619004e+03,
           1.7919020077753438063e+06,  1.7484377180034121023e+07,
           2.4086549640874004605e+04,  1.9892311916950979121e+05,
           1.8370737967633072978e+08,  2.0679040329451551508e+09,
           2.4827519375935888472e+10,  3.1669454981734887315e+11,
           4.2771126865134715582e+12,  6.0971132411392560749e+13,
           9.1486942234356396792e+14,  1.4413525170009350101e+16,
           2.3788844395175757942e+17,  4.1046081600946921885e+18,
           7.3900049415704853993e+19,  1.3859220004603943141e+21,
           2.7030825930275761623e+22,  5.4747478619645573335e+23,
           1.1498937014386333524e+25,  2.5014180692753603969e+26]

    zmag = abs(z)
    if zmag < 6
        mh1 = modhankel1(z)
        mh2 = modhankel2(z)
        mh1p = modhankel1prime(z)
        mh2p = modhankel2prime(z)
        return mh1, mh2, mh1p, mh2p
    else
        # Asymptotic expansion
        α = complex(8.53667218838951e-1)  # XXX: no idea what this is
        zpower = complex(1)
        mpower = complex(1)
        sum1 = complex(1)
        sum2 = complex(1)
        sum3 = complex(0)
        sum4 = complex(0)
        rootz = sqrt(z)
        rootz_cubed = rootz*z
        zterm = im/rootz_cubed
        mterm = -zterm
        dm = complex(0)
        term3 = complex(1)

        last = false
        m = 1
        while !last & (m <= 30)
            zpower *= zterm
            mpower *= mterm
            dm += complex(1)
            term1 = cap[m]*zpower
            term2 = cap[m]*mpower

            abs(term2/term3) >= 1 && (last = true)

            sum1 += term1
            sum2 += term2
            sum3 += term1*dm
            sum4 += term2*dm
            term4 = term2*dm

            (abs(real(term4)) <= 1e-5abs(real(sum4))) &
                (abs(imag(term4)) <= 1e-5abs(imag(sum4))) && (last = true)

            term3 = term2
            m += 1
        end
        sum3 *= zterm*complex(-1.5)/z
        sum4 *= zterm*complex(-1.5)/z

        term1 = (complex(-0.25)-im*rootz_cubed)/z
        term2 = (complex(-0.25)+im*rootz_cubed)/z

        zh1 = sum2
        zh1p = sum2*term2 + sum4
        zh2 = sum1
        zh2p = sum1*term1 + sum3

        zexp = -im*2/3*rootz_cubed + im*π*5/12

        if real(z) < 0
            exp1 = exp(zexp)
            exp2 = exp(zexp-im*π*4/3)

            if imag(z) >= 0
                zh2 *= exp1
                zh2p *= exp1

                if !bothnegative
                    zh2 += zh1/exp2
                    zh2p += zh1p/exp2
                end

                zh1 /= exp1
                zh1p /= exp1
            else
                th2 = -zh1/exp2
                th2p = -zh1p/exp2

                zh1 = zh1/exp1 + zh2*exp2
                zh1p = zh1p/exp1 + zh2p*exp2

                if !bothnegative
                    zh2 *= exp1
                    zh2p *= exp1
                else
                    zh2 = th2
                    zh2p = th2p
                end
            end
            zexp = complex(0)
        end

        zterm = α/sqrt(rootz)
        zh1 *= zterm
        zh2 *= zterm
        zh1p *= zterm
        zh2p *= zterm
    end

    return zh1, zh2, zh1p, zh2p
end


z = complex(24.6770630,2.7361517)
h1 = complex(0.3822204,-0.0108716)
h2 = complex(0.3823441,-0.0102396)
h1p = complex(-0.0548274,1.9051478)
h2p = complex(0.0503774,-1.9045298)

#  Note: this function is not currently in use
mh1test, mh2test, mh1ptest, mh2ptest = ModeFinder.modhankel(z, false)

@test mh1test ≈ h1 atol=1e-6
@test mh2test ≈ h2 atol=1e-6
@test mh1ptest ≈ h1p atol=1e-6
@test mh2ptest ≈ h2p atol=1e-6
