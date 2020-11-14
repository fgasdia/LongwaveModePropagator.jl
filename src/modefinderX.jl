#==
Functions related to the modified mode equation formulated by Morfitt and Shellman (1976) as
part of their `MODESRCH` procedure. The modified mode equation has identical zeros to the
physical mode equation except it has no poles in the complex plane
_at the reflection height_. At the ground, the poles appear identical to the physical mode
equation.

In its current form, BPM always references the ground and ionosphere reflection coefficients
at the ground height (0), so there is no advantage in using these modified equations
compared to the physical equations.
==#

# TODO: include reference height field
struct ModifiedModeEquation{W<:HomogeneousWaveguide} <: ModeEquation
    ea::EigenAngle
    frequency::Frequency
    waveguide::W
end
ModifiedModeEquation(f::Frequency, w::HomogeneousWaveguide) =
    ModifiedModeEquation(EigenAngle(complex(0.0)), f, w)
setea(ea::EigenAngle, p::ModifiedModeEquation) = ModifiedModeEquation(ea, p.frequency, p.waveguide)

"""
    R2X(ea::EigenAngle, R)

Convert reflection matrix `R` to modified reflection matrix `X` for `EigenAngle` `ea` where
``X = (R + I)/cosθ``.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
@inline function R2X(ea::EigenAngle, R)
    Cinv = ea.secθ
    return (R + I)*Cinv
end

"""
    X2R(ea::EigenAngle, X)

Convert modified reflection matrix `X` to reflection matrix `R` for `EigenAngle` `ea` where
``R = cosθ*X - I``.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
@inline function X2R(ea::EigenAngle, X)
    C = ea.cosθ
    return X2R(C, X)
end

"""
    X2R(C, X)

Convert modified reflection matrix `X` to reflection matrix `R` at cosine angle `C`.
"""
@inline X2R(C, X) = C*X - I

"""
    dXdz(X, params, z)

Return the differential of the modified reflection matrix `X` wrt height `z`.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function dXdz(X, modeequation, z)
    @unpack ea, frequency, waveguide = modeequation

    k = frequency.k
    C = ea.cosθ
    C² = ea.cos²θ

    M = susceptibility(z, frequency, waveguide)
    T = tmatrix(ea, M)

    # Precompute
    C2 = 2*C
    CT41 = C*T[4,1]

    a = SMatrix{2,2}(4*T[4,1], 0, 0, 4)
    b = SMatrix{2,2}(2*(T[4,4] - CT41), 0, -2*T[4,2], -C2)
    c = SMatrix{2,2}(2*(T[1,1]-CT41), 2*T[3,1], 0, -C2)
    d = SMatrix{2,2}(C*(T[1,1] - T[4,4]) - T[1,4] + C²*T[4,1], -C*T[3,1] + T[3,4],
            T[1,2] + C*T[4,2], C² - T[3,2])

    return -1im/2*k*(a + b*X + X*c + X*d*X)
end

function dXdCdz(XdXdC, modeequation, z)
    @unpack ea, frequency, waveguide = modeequation

    k = frequency.k
    C = ea.cosθ
    C² = ea.cos²θ

    M = susceptibility(z, frequency, waveguide)
    T = tmatrix(ea, M)
    dT = tmatrix(ea, M, DC())

    X = XdXdC[SVector(1,2),:]
    dXdC = XdXdC[SVector(3,4),:]

    # Precompute
    CT41 = C*T[4,1]

    b = SMatrix{2,2}(2*(T[4,4] - CT41), 0, -2*T[4,2], -C2)
    c = SMatrix{2,2}(2*(T[1,1]-CT41), 2*T[3,1], 0, -C2)
    d = SMatrix{2,2}(C*(T[1,1] - T[4,4]) - T[1,4] + C²*T[4,1], -C*T[3,1] + T[3,4],
            T[1,2] + C*T[4,2], C² - T[3,2])

    db = SMatrix{2,2}(2*dT[4,4], 0, 0, -2)
    dc = SMatrix{2,2}(2*dT[1,1], 0, 0, -2)
    dd = SMatrix{2,2}(C*(dT[1,1] - dT[4,4]) + T[1,1] - T[4,4] + 2*CT41 - dT[1,4],
            dT[3,4] - T[3,1], dT[1,2] + T[4,2], 2*C - dT[3,2])

    return -1im/2*k*(db*X + b*dX + dX*c + X*dc + dX*d*X + X*dd*X + X*d*dX)
end

function integratedreflectionX(modeequation::ModifiedModeEquation,
    params::IntegrationParams=DEFAULT_INTEGRATIONPARAMS) where T

    @unpack tolerance, solver, force_dtmin = params

    Mtop = susceptibility(TOPHEIGHT, modeequation)
    Rtop = sharpboundaryreflection(modeequation.ea, Mtop)
    Xtop = R2X(modeequation.ea, Rtop)

    prob = ODEProblem{false}(dXdz, Xtop, (TOPHEIGHT, BOTTOMHEIGHT), modeequation)

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, solver, abstol=tolerance, reltol=tolerance,
                force_dtmin=force_dtmin, dt=1,
                save_on=false, save_start=false, save_end=true)

    X = sol[end]

    return X
end

"""
    fresnelreflectionX(ea::EigenAngle, ground::Ground, frequency::Frequency)

Return modified ground reflection matrices `n` and `d`, corresponding to the modified
ionosphere reflection matrix `X`.

The coefficients `n` and `d` are defined at the ground and satisfy the relationship
``nd⁻¹ = (Rg⁻¹ + I)/cosθ``.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function fresnelreflectionX(ea::EigenAngle, ground::Ground, frequency::Frequency)
    C, S² = ea.cosθ, ea.sin²θ
    ω = frequency.ω

    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*E0))

    CNg² = C*Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)
    invsqrtNg²mS² = inv(sqrtNg²mS²)

    n11 = 1
    n22 = invsqrtNg²mS²

    d11 = (C - sqrtNg²mS²/Ng²)/2
    d22 = (C*invsqrtNg²mS² - 1)/2

    n = SDiagonal(n11, n22)
    d = SDiagonal(d11, d22)

    return n, d
end
fresnelreflection(m::ModifiedModeEquation) =
    fresnelreflectionX(m.ea, m.waveguide.ground, m.frequency)

function fresnelreflectionX(z, ea::EigenAngle, ground::Ground, frequency::Frequency)
    n, d = fresnelreflectionX(ea, ground, frequency)

    # If altitude `z` is ≈ 0, don't bother integrating upward

    if z > 1
        n, d = fresnelintegration(z, n, d, ea, frequency)
    end

    return n, d
end
fresnelreflection(z, m::ModifiedModeEquation) =
    fresnelreflectionX(z, m.ea, m.waveguide.ground, m.frequency)

function fresnelintegration(z, nb, db, ea::EigenAngle, frequency::Frequency)
    C, C² = ea.cosθ, ea.cos²θ

    α = 2/EARTHRADIUS
    K = 1im*cbrt(α/frequency.k)
    L = 1im*(α/(2*K))
    koatt = (frequency.k/α)^(2/3)

    n² = 1 + α*(z - CURVATURE_HEIGHT)
    nb² = 1 - α*CURVATURE_HEIGHT

    ab = C + L/nb²
    a = C + L/n²
    kb = K/nb²
    k = K/n²

    q = koatt*(C² + α*(z - CURVATURE_HEIGHT))
    qb = koatt*(C² - α*CURVATURE_HEIGHT)

    h₁, h₂, h₁p, h₂p = modifiedhankel(q)
    h₁b, h₂b, h₁pb, h₂pb = modifiedhankel(qb)

    f = h₁b*h₂ - h₂b*h₁
    fb = h₁pb*h₂ - h₂pb*h₁
    fz = h₁b*h₂p - h₂b*h₁p
    fbz = h₁pb*h₂p - h₂pb*h₁p

    d11 = -nb[1,1]*(ab*a*f + k*ab*fz + kb*a*fb + kb*k*fbz) + 2*db[1,1]*(a*f + k*fz)
    d22 = -nb[2,2]*(C²*f + K*C*fz + K*C*fb + K^2*fbz) + 2*db[2,2]*(C*f + K*fz)
    n11 = -2*nb[1,1]*(ab*f + kb*fb) + 4*db[1,1]*f
    n22 = -2*nb[2,2]*(C*f + K*fb) + 4*db[2,2]*f

    n = SDiagonal(n11, n22)
    d = SDiagonal(d11, d22)

    return n, d
end

"""
    modalequationX(X, n, d)

Return value of the modified mode equation [^Morfitt1976] given reflection coefficient `X`
for the ionosphere and ground coefficients `n` and `d`.

The modified variables ``X = (R + 1)/cosθ`` and ``nd⁻¹ = (Rg⁻¹ + 1)/cosθ``. The modified
mode equation is ``F₃(θ) = (∥n∥ - ∥X∥ ∥d∥)(⟂n⟂ - ⟂X⟂ ⟂d⟂) - ∥X⟂ ⟂X∥ ∥d⟂ ⟂d∥``. The modified
equation is not the same as the original, but except for the point where θ = 90°, they have
the same zeros.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function modalequationX(X, n, d)
    return (n[1,1] - X[1,1]*d[1,1])*(n[2,2] - X[2,2]*d[2,2]) - X[1,2]*X[2,1]*d[1,1]*d[2,2]
end

function solvemodalequation(z, modeequation::ModifiedModeEquation,
    integrationparams::IntegrationParams=DEFAULT_INTEGRATIONPARAMS)

    X = integratedreflectionX(modeequation, integrationparams)
    n, d = fresnelreflection(z, modeequation)

    f = modalequationX(X, n, d)
    return f
end
function solvemodalequation(θ, z, modeequation::ModifiedModeEquation,
    integrationparams::IntegrationParams=DEFAULT_INTEGRATIONPARAMS)
    # Convenience function for `grpf`
    modeequation = setea(EigenAngle(θ), modeequation)
    solvemodalequation(z, modeequation, integrationparams)
end

"""
LWPC version of MS76 pg 23-24
"""
function freespaceintegration(z, Xb, ea::EigenAngle, frequency::Frequency; zb=BOTTOMHEIGHT)
    C, C² = ea.cosθ, ea.cos²θ

    α = 2/EARTHRADIUS
    K = 1im*cbrt(α/frequency.k)
    L = 1im*(α/(2*K))
    koatt = (frequency.k/α)^(2/3)

    n² = 1 + α*(z - CURVATURE_HEIGHT)
    nb² = 1 + α*(zb - CURVATURE_HEIGHT)

    ab = C + L/nb²
    a = C + L/n²
    kb = K/nb²
    k = K/n²

    q = koatt*(C² + α*(z - CURVATURE_HEIGHT))
    qb = koatt*(C² + α*(zb - CURVATURE_HEIGHT))

    h₁, h₂, h₁p, h₂p = modifiedhankel(q)
    h₁b, h₂b, h₁pb, h₂pb = modifiedhankel(qb)

    f = h₁b*h₂ - h₂b*h₁
    fb = h₁pb*h₂ - h₂pb*h₁
    fz = h₁b*h₂p - h₂b*h₁p
    fbz = h₁pb*h₂p - h₂pb*h₁p

    ep1 = -Xb[1,1]/2*(ab*a*f + k*ab*fz + kb*a*fb + kb*k*fbz) + a*f + k*fz
    ep2 = -Xb[1,2]/2*(ab*a*f + k*ab*fz + kb*a*fb + kb*k*fbz)
    ey1 = -Xb[2,1]/2*(C²*f + K*C*fz + K*C*fb + K^2*fbz)
    ey2 = -Xb[2,2]/2*(C²*f + K*C*fz + K*C*fb + K^2*fbz) + C*f + K*fz
    den = ep1*ey2 - ep2*ey1

    v1 = -Xb[1,1]*(ab*f + kb*fb) + 2*f
    v2 = -Xb[1,2]*(ab*f + kb*fb)
    v3 = -Xb[2,1]*(C*f + K*fb)
    v4 = -Xb[2,2]*(C*f + K*fb) + 2*f

    X11 = (v1*ey2 - v2*ey1)/den
    X12 = (v2*ep1 - v1*ep2)/den*nb²
    X21 = (v3*ey2 - v4*ey1)/den/nb²
    X22 = (v4*ep1 - v3*ep2)/den

    return SMatrix{2,2}(X11, X21, X12, X22)
end

function freespaceintegration(z, Xb, dXb, ea::EigenAngle, frequency::Frequency, ::DC;
    zb=BOTTOMHEIGHT)

    C, C² = ea.cosθ, ea.cos²θ

    α = 2/EARTHRADIUS
    K = 1im*cbrt(α/frequency.k)
    L = 1im*(α/(2*K))
    koatt = (frequency.k/α)^(2/3)

    n² = 1 + α*(z - CURVATURE_HEIGHT)
    nb² = 1 + α*(zb - CURVATURE_HEIGHT)

    ab = C + L/nb²
    a = C + L/n²
    kb = K/nb²
    k = K/n²

    q = koatt*(C² + α*(z - CURVATURE_HEIGHT))
    qb = koatt*(C² + α*(zb - CURVATURE_HEIGHT))

    dq = 2*C*koatt

    h₁, h₂, h₁p, h₂p = modifiedhankel(q)
    h₁b, h₂b, h₁pb, h₂pb = modifiedhankel(qb)

    f = h₁b*h₂ - h₂b*h₁
    fb = h₁pb*h₂ - h₂pb*h₁
    fz = h₁b*h₂p - h₂b*h₁p
    fbz = h₁pb*h₂p - h₂pb*h₁p

    df = (fb + fz)*dq
    dfb = (fbz - qb*f)*dq
    dfz = (fbz - q*f)*dq
    dfbz = -(q*fb + qb*fz)*dq

    ep1 = -Xb[1,1]/2*(ab*a*f + k*ab*fz + kb*a*fb + kb*k*fbz) + a*f + k*fz
    ep2 = -Xb[1,2]/2*(ab*a*f + k*ab*fz + kb*a*fb + kb*k*fbz)
    ey1 = -Xb[2,1]/2*(C²*f + K*C*fz + K*C*fb + K^2*fbz)
    ey2 = -Xb[2,2]/2*(C²*f + K*C*fz + K*C*fb + K^2*fbz) + C*f + K*fz
    den = ep1*ey2 - ep2*ey1

    dep1 = -dXb[1,1]/2*(ab*a*f + k*ab*fz + kb*a*fb + kb*k*fbz) -
            Xb[1,1]/2*(ab*f + a*f + ab*a*df + k*(fz + ab*dfz) + kb*(fb + a*dfb) + kb*k*dfbz) +
                f + a*df + k*dfz
    dep2 = -dXb[1,2]/2*(ab*a*f + k*ab*fz + kb*a*fb + kb*k*fbz) -
            Xb[1,2]/2*(ab*f + a*f + ab*a*df + k*(fz + ab*dfz) + kb*(fb + a*dfb) + kb*k*dfbz)
    dey1 = -dXb[2,1]/2*(C²*f + K*C*fz + K*C*fb + K^2*fbz) -
            Xb[2,1]/2*(2*C*f + C²*df + K*(fz + C*dfz) + K*(fb + C*dfb) + K^2*dfbz)
    dey2 = -dXb[2,2]/2*(C²*f + K*C*fz + K*C*fb + K^2*fbz) + f + C*df + K*dfz -
            Xb[2,2]/2*(2*C*f + C²*df + K*(fz + C*dfz) + K*(fb + C*dfb) + K^2*dfbz)
    dden = ep1*dey2 + dep1*ey2 - ep2*dey1 - dep2*ey1

    v1 = -Xb[1,1]*(ab*f + kb*fb) + 2*f
    v2 = -Xb[1,2]*(ab*f + kb*fb)
    v3 = -Xb[2,1]*(C*f + K*fb)
    v4 = -Xb[2,2]*(C*f + K*fb) + 2*f

    dv1 = -dXb[1,1]*(ab*f + kb*fb) - Xb[1,1]*(f + ab*df + kb*dfb) + 2*df
    dv2 = -dXb[1,2]*(ab*f + kb*fb) - Xb[1,2]*(f + ab*df + kb*dfb)
    dv3 = -dXb[2,1]*(C*f + K*fb) - Xb[2,1]*(f + C*df + K*dfb)
    dv4 = -dXb[2,2]*(C*f + K*fb) - Xb[2,2]*(f + C*df + K*dfb) + 2*df

    X11 = (v1*ey2 - v2*ey1)/den
    X12 = (v2*ep1 - v1*ep2)/den*nb²
    X21 = (v3*ey2 - v4*ey1)/den/nb²
    X22 = (v4*ep1 - v3*ep2)/den

    dX11 = (v1*dey2 + dv1*ey2 - v2*dey1 - dv2*ey1)/den - Xb[1,1]*dden/den
    dX12 = (v2*dep1 + dv2*ep1 - v1*dep2 - dv1*ep2)/den*nb² - Xb[1,2]*dden/den
    dX21 = (v3*dey2 + dv3*ey2 - v4*dey1 - dv4*ey1)/den/nb² - Xb[2,1]*dden/den
    dX22 = (v4*dep1 + dv4*ep1 - v3*dep2 - dv3*ep2)/den - Xb[2,2]*dden/den

    return SMatrix{2,2}(X11, X21, X12, X22), SMatrix{2,2}(dX11, dX21, dX12, dX22)
end

"""
Described in MS76 pg 16.

in LWPC, this is mostly mf_initlg.for or mf_initev.for if that fails
"""
function reflectionheight(modeequation::ModifiedModeEquation,
    params::IntegrationParams=DEFAULT_INTEGRATIONPARAMS) where T

    @unpack tolerance, solver, force_dtmin = params
    frequency = modeequation.frequency

    altitudes = 40e3:1e3:100e3
    eas = SVector(EigenAngle(deg2rad(70.0-0.1im)), EigenAngle(deg2rad(70.0-8im)),
                  EigenAngle(deg2rad(89.9-0.1im)), EigenAngle(deg2rad(89.9-8im)))

    Xbs = MMatrix{length(eas),SMatrix{2,2,ComplexF64,4}}(undef)
    for t in eachindex(eas)
        modeequation = setea(eas[t], modeequation)
        Xbs[t] = integratedreflectionX(modeequation)
    end

    minXsum = Inf
    minXaltitude = Inf
    for i in eachindex(altitudes)
        Xsum = 0.0
        for t in eachindex(eas)
            X = freespaceintegration(altitudes[i], Xbs[t], eas[t], frequency)
            Xsum += abs2(X[1,1]) + abs2(X[2,1]) + abs2(X[1,2]) + abs2(X[2,2])
        end

        if Xsum < minXsum
            minXsum = Xsum
            minXaltitude = altitudes[i]
        end
    end

    return minXsum, minXaltitude
end
