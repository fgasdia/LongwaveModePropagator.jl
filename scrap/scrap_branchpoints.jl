
"""
TEMP
"""
function branchpoints(frequency::Frequency, bfield::BField, species::Species)
    B, lx, ly, lz = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    lx², ly², lz² = lx^2, ly^2, lz^2
    ω = frequency.ω

    # Constitutive relations (see Budden1955a, pg. 517 or Budden1988 pg. 39)
    e, m, N, ν = species.charge, species.mass, species.numberdensity, species.collisionfrequency
    invω = inv(ω)
    invmω = invω/m  # == inv(m*ω)

    z = TOPHEIGHT
    Y = abs(e*B*invmω)  # |ωₕ/ω|  gyrofrequency / ω  # Nagano et al 1975 specifies |ωₕ/ω|
    Z = ν(z)*invω  # collision frequency / ω

    g = Y/(Z + 1im)

    # maybe?
    ec = 2*(50e3 - TOPHEIGHT)/EARTHRADIUS
    p² = 1 + ec
    p = sqrt(p²)

    t = [p*(1 + sqrt(1 + g^2))/g, p*(1 - sqrt(1 + g^2))/g]

    bpt = Vector{ComplexF64}(undef, 6)
    qq = similar(bpt)

    i = 1
    for ii = 1:2
        rt = lx*sqrt((lx² + lz²)*p² - t[ii]^2)
        qp = (lz*t[ii] + rt)/(lx² + lz²)

        for nrt = 1:2
            q² = qp^2
            s = sqrt(p² - q²)
            real(s) < 0 && (s = -s)

            a = lx*s + lz*qp
            b = -lz*s + lx*qp
            z² = ((b^2 + ly²*p²)*g)^2 - 4*a^2*p²
            z²mgp = abs(z²)

            qm = -qp
            a = lx*s + lz*qm
            b = -lz*s + lx*qm
            z² = ((b^2 + ly²*p²)*g)^2 - 4*a^2*p²
            z²mgm = abs(z²)

            if z²mgp < z²mgm
                qq[i] = qp
                z²mag = z²mgp
            else
                qq[i] = qm
                z²mag = z²mgm
            end

            c² = qq[i]^2 - ec
            c = sqrt(c²)
            real(qq[i]) < 0 && (c = -c)

            tcompl = -1im*log(s + c*1im)
            bpt[i] = 90 - tcompl*57.296  # NOTE IN DEGREES!

            i += 1
            qp = (lz*t[ii] - rt)/(lx² + lz²)
        end
    end

    bpt[5] = complex(90)

    if ec < 0
        qpt = complex(acos(sqrt(-ec))*57.296)
    else
        argln = sqrt(ec) + sqrt(1 + ec)
        qpt = 90 + -1im*log(argln)*57.296
    end
    bpt[6] = qpt

    return bpt
end
