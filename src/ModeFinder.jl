module ModeFinder

thisdir = dirname(@__FILE__())
any(path -> path==thisdir, LOAD_PATH) || push!(LOAD_PATH, thisdir)

using SpecialFunctions
using LinearAlgebra
using StaticArrays
using PolynomialRoots
using DifferentialEquations

using LWMS
using IonosphereProfile


const speedoflight = 299792458  # m/s
const μ₀ = 4e-7π  # H/m
const ϵ₀ = 1/(μ₀*speedoflight^2)  # F/m
const mₑ = 9.10938356e-31  # kg
const fundamentalcharge = 1.602176621e-19  # C


mutable struct Mode
    theta::ComplexF64
    c::ComplexF64
    s::ComplexF64
    csq::ComplexF64
    ssq::ComplexF64
    omega::Float64
    wavenumber::Float64
    ideriv::Bool
end

mutable struct Area
    gmax::Float64
    boxratio::Float64
    meshangle::Float64
    autosubdivide::Bool
end

mutable struct ConstitutiveRelations
    X::Float64
    Y::Float64
    Z::Float64
    U::ComplexF64
end

mutable struct DirectionCosines
    dcl::Float64
    dcm::Float64
    dcn::Float64
    G::ComplexF64
end


"""
MF_DRIVER sets omega, wavenr, gmax, tmesh, lub, ranger, rangei, and calls MF_WVGD.

LWPC's MF_DRIVER is a short subroutine called by LWP_DRIVER to get starting solutions for
modesearch in an individual path segment. I believe _starting solutions_ are the approximate
eigenangles. If eigenangles (_modes_) are found, then LWP_DRIVER proceeds to enter SW_STEP
and additional routines, presumably to find the exact or complete solutions.
"""
function driver(inputs::Inputs, modes::LWMSModes)
# Unpack
freq = inputs.freq

ω = 2π*freq
k = ω/speedoflight
meshangle = sqrt(1/freq)
modes.deltathetathreshold = sqrt(15/freq)/100

if freq <= 29.9e3
    gmax = 3.0
elseif freq <= 99.9e3
    gmax = 2.0
else
    gmax = 1.0
end

if modes.ranger[1] > modes.ranger[2]
    modes.ranger[1], modes.ranger[2] = modes.ranger[2], modes.ranger[1]
end
if modes.rangei[1] < modes.rangei[2]
    modes.rangei[1], modes.rangei[2] = modes.rangei[2], modes.rangei[1]
end

mode = Mode(complex(NaN, NaN),
            complex(NaN, NaN), complex(NaN, NaN),
            complex(NaN, NaN), complex(NaN, NaN),
            ω, k, false)
area = Area(gmax, NaN, meshangle, false)

waveguide(mode)

return
end

"""
LWPC's MF_WVGD is a medium length subroutine that obtains approximate eigenangles (_modes_)
in a single path segment. It calls several secondary routines to accomplish this task.
"""
function waveguide(inputs::Inputs)
# Unpack
azim, dip = inputs.azim, inputs.dip

# Initialize
dcl = cosd(dip)*cosd(azim)
dcm = cosd(dip)*sind(azim)
dcn = -sind(dip)

initializeeigens()
boundarys()

# >>>
end

"""
MF_INITEV initializes the eigenangle search by calculating the profile, some basic
magnetoionic parameters, and calling MF_BRNCPT if necessary. It also determines htntrp.

LWPC's MF_INITEV is a short subroutine. It determines `htntrp`, which is sometimes used as a
reference height. It evaluates the profile and calcultes `x`, `y`, `z`, and `g`. If `g` is
not low, it calls MF_BRNCPT to locate any branch points in the search grid.

The height HTNTRP at which cap X (as defined by Budden) is equal to the input or default
value of X = Xinterp. It is used as a default value of the reference height, D, in LAGRNG in
those cases in which the usual procedure of that routine for finding an appropriate
reference height cannot be used and REFHT is not used. It is always the reference height
used in MF_INTROE.

From NOSC TR1143:
- `x = cx*en[1]` is the normalized electron density
- `y = -cy*bfield/omega` is the normalized magnetic field strength
- `z = nu[1]/omega` is ν/ω (collision frequency/angular wave frequency)
- `g` is ``y/(z+i)`` (the ratio of magnetic field to collision frequency)

Some additional values used in MF_INITEV that are not obvious:
- `cx` is in Morfitt Shellman 1976 pg 77 eq 11, which references Wait's NBS TN 300 where
```math
ωᵣ(z) = 2.5\times 10^5 exp(β(z-h′)) = 3.182357\times 10^9 N(z)/ν(z) \\
ω₀(z) = 3.18\times 10^9 N(0)
```
and `cx=3.182357e9`, ``ω₀`` is the plasma frequency, ``N(z)`` is the electron density in electrons/cubic centimeter, and
``ν(z)`` is the collision frequency in collisions/second.
- `cy = 1.758796e7` is derived (by me) from the ``Y`` in Morfitt Shellman 1976 pg 14 where
``Y = μ₀ ωₘ/ω`` henries where ωₘ is the magnetic gyrofrequency (1/s)
so that `cy` here is
```math
-cy*bfield/ω = μ₀ ωₘ/ω = μ₀ (eB/mc)/ω
cy = μ₀ (e/mc) * (B/ω)
```
except that 1.758796e7 is just e/mc, so it appears to have no μ₀.
"""
function initialheightinterp(ω, height, topheight, bottomheight, xinterp,
                             spec::Constituent)
    e, m, N, ν = spec.charge, spec.mass, spec.numberdensity, spec.collisionfrequency
    X = N(height)*e^2/(ϵ₀*m*ω^2)
    Y = e*B/(m*ω)
    Z = ν(height)/ω
    U = 1 - im*Z

    ht = bottomheight
    x = 0.0
    while x < xinterp
        prfl_ennu(ht, en, nu)  # TODO: Can we feed this arrays to avoid the loop?
        x = coeffen*en[1]
        if ht == topheight
            error("`htinterp` cannot be set")
        end
        ht += 1
    end

    heights.htinterp = ht - 2.0
    prfl_ennu(htinterp, en, nu)

    z = nu[1]/ω
    y = -cy*bfield/ω
    drcs.g = y/(z + im)
    drcs.ec = 2(heights.htinterp - h)/earthradius

    # If the magnitude of G is very small, that is, if the collision frequency at HTNTRP is
    # large, MF_INTROE is not used, that is, reflection coefficients are not used. In this
    # case is done entirely in MF_LAGRNG.
    if abs(drcs.g) < parameters.gthreshold
        flags.lowg = true
    else
        flags.lowg = false
        findbranchpoints()
    end
    return
end

"""
Computation of susceptibility `M` and temporary `T` matrices as defined by Budden (1955)
[](10.1098/rspa.1955.0027).

Constitutive relations (see Budden 1955, pg. 517):
```math
X = Ne²/ϵ₀mω² = ωₚ²/ω²
Y = eB/mω = ωₘ/ω
Z = ν/ω
U = 1 - iZ
```

Susceptibility matrix ``M``:
```math
\\mathbf{M} = -\\frac{X}{U(U^2-Y^2)}
    \\begin{pmatrix}
        U^2 - l^2Y^2 & -iUnY - lmY^2 & iUmY-lnY^2
        iUnY - lmY^2 & U^2 -m^2Y^2 & -iUlY-mnY^2
        -iUmY-lnY^2 & iUlY-mnY^2 & U^2-n^2Y^2
    \\end{pmatrix}
```

Matrix ``T`` fulfills the differential equation ``e′ = -iTe`` where ``e`` is the column
matrix ``e = [Eₓ -Ey Hₓ Hy]``.

Note:
- `T44` in MS 76 should be negative
"""
function tmatrix(θ, height, referenceheight, ω, spec::Constituent, B, dcl, dcm, dcn)
    # Initialize
    S = sind(θ)
    C = cosd(θ)
    earthcurvature = 2/earthradius*(height - referenceheight)

    # TODO: Add M up for each species
    # Constitutive relations (see Budden 1955, pg. 517)
    e, m, N, ν = spec.charge, spec.mass, spec.numberdensity, spec.collisionfrequency
    X = N(height)*e^2/(ϵ₀*m*ω^2)
    Y = e*B/(m*ω)
    Z = ν(height)/ω
    U = 1 - im*Z

    @debug "Constitutive relations: " X Y Z U

    # Construct susceptibility matrix (see Budden 1955, eq. 3)
    M11 = U^2 - (dcl*Y)^2
    M22 = U^2 - (dcm*Y)^2
    M33 = U^2 - (dcn*Y)^2

    M12 = -im*dcn*Y*U - dcl*dcm*Y^2
    M21 = im*dcn*Y*U - dcl*dcm*Y^2

    M13 = im*dcm*Y*U - dcl*dcn*Y^2
    M31 = -im*dcm*Y*U - dcl*dcn*Y^2

    M23 = -im*dcl*Y*U - dcm*dcn*Y^2
    M32 = im*dcl*Y*U - dcm*dcn*Y^2

    M = @MMatrix [M11 M12 M13;
                  M21 M22 M23;
                  M31 M32 M33]

    M .*= -X/(U*(U^2 - Y^2))
    # XXX: In LWPC `earthcurvature` is not multiplied by capd (the above line)
    # Is ^this^ correct? It _is_ multiplied in MS 1976
    M[1,1] += earthcurvature
    M[2,2] += earthcurvature
    M[3,3] += earthcurvature

    oneplusM33 = 1 + M[3,3]
    T = @SMatrix [-S*M[3,1]/oneplusM33 S*M[3,2]/oneplusM33 0 (C^2+M[3,3])/oneplusM33;
                  0 0 1 0;
                  M[2,3]*M[3,1]/oneplusM33-M[2,1] C^2+M[2,2]-M[2,3]*M[3,2]/oneplusM33 0 S*M[2,3]/oneplusM33;
                  1+M[1,1]-M[1,3]*M[3,1]/oneplusM33 M[3,2]*M[1,3]/oneplusM33-M[1,2] 0 -S*M[1,3]/oneplusM33]

    return M, T
end

"""
Compute the (S matrix) coefficients used in the differential equations for ``(R+1)/C``.

The coefficients are derived from the differential equation for the reflection
coefficient matrix given by Budden (1955) [](10.1098/rspa.1955.0027). This function returns
`A`, `B`, `C`, and `D` which are a rearranged form of the S matrix components for X, rather
than R (see Morfitt Shellman 1976).

The differential equations are of the form
```math
dR/dz = -ik/2 \\left( S^{(21)} + S^{(22)}R + RS^{(11)} - RS^{(12)}R \\right)
```
where `S` is constructed from the `T`, and therefore `M` (susceptibility), matrix. It was
rearranged to form the differential equations in terms of ``X = (R+1)/C``, giving
```math
dX/dz = -ik/2 \\left( A + BX + XC + XDX \\right)
```
where ``k`` is the wave number, ``z`` is the height, and
```math
\\begin{align}
A &= \\begin{pmatrix}
        4T_{41} & 0 \\
        0 & 4
    \\end{pmatrix} \\
B &= 2\\begin{pmatrix}
        T_{44} & -T_{42} \\
        0 & -C
    \\end{pmatrix} \\
C &= 2\\begin{pmatrix}
        T_{11}-T_{44} & 0 \\
        T_{31} & -C
    \\end{pmatrix} \\
D &= \\begin{pmatrix}
        CT_{11}-T_{44}-T_{14}+C^2T_{41} & T_{12}+CT_{42} \\
        -CT_{31}+T_{34} & C^2-T_{32}
    \\end{pmatrix}
\\end{align}
```
Note:
- MS 76' says ``C = -S^{11} + S^{12}`` and then writes that in terms of ``T``, but `C11` is
missing a minus sign, it should be ``-T^{11} - CT^{41}``.
"""
function smatrix(θ, T)
    c = cosd(θ)

    A = @SMatrix [4T[4,1]   0;
                  0         4]
    B = @SMatrix [2(T[4,4]-c*T[4,1])    -2T[4,2];
                  0                     -2c]
    C = @SMatrix [-2(T[1,1]+c*T[4,1])   0;
                  2T[3,1]               -2c]
    D = @SMatrix [c*(T[1,1]-T[4,4])-T[1,4]+c^2*T[4,1]   T[1,2]+c*T[4,2];
                  -c*T[3,1]+T[3,4]                      c^2-T[3,2]]

    return A, B, C, D
end

"""
Computation of ``X = (R+1)/C`` for reflection from a sharply bounded anisotropic ionosphere of
semi-infinite extent where R is the reflection coefficient matrix as defined by Budden. The
solution is used as the initial condition for integration. This solution is derived from the
one in Radio Science Vol. 3, Aug 1968 pg. 792-795, although makes use of a different solver.

The physics is described by ``D\\mathbf{E} = 0`` where ``E`` is the electric vector of the
elm wave and ``D`` is the dispersion matrix ``I + M + L`` where ``I`` is the unit matrix and
```math
L = \\begin{pmatrix}
        -q^2 & 0 & qS \\
        0 & -q^2 - S^2 & 0 \\
        qS & 0 & -S^2
    \\end{pmatrix}
q = nₜ\\cos θₜ
θₜ = \\mathrm{angle of refraction}
S = \\sin θ
θ = \\mathrm{angle of incidence}
```
If ``D\\mathbf{E} = 0`` is to have nontrivial solutions, the determinant of ``D`` must be zero.
The equation ``|D| = 0`` is a 4th order polynomial ("quartic") in `q`. The two roots of the
quartic closest to the positive real and negative imaginary axis are calculated and chosen
to form ``D``, which are then used to form the reflection coefficient matrix ``R`` (or ``X``).
"""
function sharplyboundedreflectionmatrix(θ, M)
    # Initialize
    S = sind(θ)
    C = cosd(θ)
    P = MVector{2,ComplexF64}(undef)
    T = MVector{2,ComplexF64}(undef)

    # Booker quartic coefficients
    b4 = 1 + M[3,3]
    b3 = S*(M[1,3] + M[3,1])
    b2 = -(C^2 + M[3,3])*(1 + M[1,1]) + M[1,3]*M[3,1] -
        (1 + M[3,3])*(C^2 + M[2,2]) + M[2,3]*M[3,2]
    b1 = S*(M[1,2]*M[2,3] + M[2,1]*M[3,2] -
        (C^2 + M[2,2])*(M[1,3] + M[3,1]))
    b0 = (1 + M[1,1])*(C^2 + M[2,2])*(C^2 + M[3,3]) +
        M[1,2]*M[2,3]*M[3,1] + M[1,3]*M[2,1]*M[3,2] -
        M[1,3]*(C^2 + M[2,2])*M[3,1] -
        (1 + M[1,1])*M[2,3]*M[3,2] -
        M[1,2]*M[2,1]*(C^2 + M[3,3])

    q = roots([b0, b1, b2, b3, b4])

    # Calculate angle from 315°
    function anglefrom315(qval)
        angq = rad2deg(angle(qval))
        angq < 0 && (angq += 360)
        angq < 135 && (angq += 360)
        abs(angq - 315)
    end
    sort!(q, by=anglefrom315)

    # Dispersion matrix elements for R, X
    D12 = M[1,2]
    D32 = M[3,2]
    D33 = C^2 + M[3,3]
    for j = 1:2
        D11 = 1 + M[1,1] - q[j]^2
        D13 = M[1,3] + q[j]*S
        D31 = M[3,1] + q[j]*S
        denom = D11*D33 - D13*D31
        P[j] = (-D12*D33 + D13*D32)/denom
        T[j] = q[j]*P[j] - S*(-D11*D32 + D12*D31)/denom
    end
    Δ = (T[1]*C + P[1])*(C + q[2]) - (T[2]*C + P[2])*(C + q[1])

    # XXX: Do ∥X∥ and ⟂X⟂ really equal (∥R∥ + 1)/C and (⟂R⟂ + 1)/C? Using `r11` and `r22` from
    # mf_initr.for for now, but I haven't been able to prove this is true mathematically.
    X = SMatrix{2, 2}(2(T[1]*(C + q[2]) - T[2]*(C +q[1]))/Δ,
                      -2(q[1] - q[2])/Δ,
                      -2(T[1]*P[2] - T[2]*P[1])/Δ,
                      2((T[1]*C + P[1]) - (T[2]*C + P[2]))/Δ)

    return q, X
end

"""
Return the differential equation that describes the modified reflection coefficient as a
function of height, ``dX/dz``.

The differential equations are of the form
```math
dR/dz = -ik/2 \\left( S^{(21)} + S^{(22)}R + RS^{(11)} - RS^{(12)}R \\right)
```
where `S` is constructed from the `T`, and therefore `M` (susceptibility), matrix. It was
rearranged to form the differential equations in terms of ``X = (R+1)/C``, giving
```math
dX/dz = -ik/2 \\left( A + BX + XC + XDX \\right)
```
"""
function dXdh(X, A, B, C, D, k)
    # TODO: An in place version of this?
    (-im*k/2)*(A + B*X + X*C + X*D*X)
end

"""
"""
function xderivative(X, params, height)
    T = tmatrix(params.θ, height, params.referenceheight, params.ω,
                params.species, params.B, params.dcl, params.dcm, params.dcn)[2]
    A, B, C, D = smatrix(params.θ, T)
    dXdh(X, A, B, C, D, params.k)
end

"""
Full wave solution for `X` through ionosphere.

Integrates ``dX/dh`` (`dXdh()`) from `topheight` to `bottomheight`.
"""
function integratethroughionosphere(θ, ω, topheight, bottomheight, referenceheight,
                                    spec, B, dcl, dcm, dcn)
    # Initial condition for integration
    M = tmatrix(θ, topheight, referenceheight, ω, spec, B, dcl, dcm, dcn)[1]
    initialX = sharplyboundedreflectionmatrix(θ, M)[2]

    k = ω/speedoflight*1e3  # k must be based in km because we integrate over height in km
    params = (θ=θ, referenceheight=referenceheight, ω=ω, k=k, species=spec,
              B=B, dcl=dcl, dcm=dcm, dcn=dcn)

    heightbounds = (topheight, bottomheight)
    prob = ODEProblem(xderivative, initialX, heightbounds, params)
    sol = solve(prob, BS3(), reltol=1e-3)  # Faster, less accurate
    # sol = solve(prob, Tsit5(), reltol=1e-4)  # Decent accuracy for the speed

    return sol
end

"""
These functions calculate the modified Hankel functions of order 1/3 and their derivatives.

The functions h₁(ζ) and h₂(ζ) satisfy the Stokes differential equation (_Airy_ function)
```math
\\frac{d²w}{dζ²} + ζw = 0
```
"""
modhankel1(z) = 12^(1/6)*exp(-π/6*im)*(SpecialFunctions.airyai(-z)-im*SpecialFunctions.airybi(-z))
modhankel2(z) = 12^(1/6)*exp(π/6*im)*(SpecialFunctions.airyai(-z)+im*SpecialFunctions.airybi(-z))
modhankel1prime(z) = Complex(-1)^(5/6)*12^(1/6)*(SpecialFunctions.airyaiprime(-z)-im*SpecialFunctions.airybiprime(-z))
modhankel2prime(z) = Complex(-1)^(-5/6)*12^(1/6)*(SpecialFunctions.airyaiprime(-z)+im*SpecialFunctions.airybiprime(-z))

"""
This function calculates the modified Hankel functions of order 1/3 and their derivatives.

The functions h₁(ζ) and h₂(ζ) satisfy the Stokes differential equation (_Airy_ function)
```math
\\frac{d²w}{dζ²} + ζw = 0
```

Note: This may actually be doing more than just calculating the hankel functions, because
when used with `integratethroughfreespace()` below, which follows the math steps in MS 1976,
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

"""
Performs an integration of the differential equations for the ionosphere reflection matrix
through a free space region over a curved earth.

The integration may be performed in either a positive or negative height direction, but LWPC
always integrates upward. The integration variables are the matrix ``(R+1)/C`` where ``R`` is
the reflection matrix described by Budden. This solution is based on Budden, Radio Waves in
the Ionosphere, in particular the material on pg 118, 327-329, 336-338, and 343-345. MS
1976, especially Appendix 2, presents a detailed explanation.

Computation at θ = 90° is not excluded. Also, the equations are formulated so that a smooth
transition is made from the "curved earth" form to the "flat earth" form for appropriate
values of θ.

From MS 1976 pg. 17 or Budden, Radio Waves in the Ionosphere, pg. 118, the total field
components just below the free space-ionosphere boundary are
```math
Eₓ = (E∥ⁱ - E∥ʳ)\\cos θ_I
E_y = E_yⁱ + E_yʳ
Hₓ = (E_yʳ - E_yⁱ)\\cos θ_I
H_y = E∥ⁱ + E∥ʳ
```
where ``θ_I`` is the angle the incident wave normal makes with the vertical ``z`` axis. Also
we note
```math
⟂R⟂ = E_yʳ / E_yⁱ
∥R⟂ = E_yʳ / E∥ⁱ
∥R∥ = E∥ʳ / E∥ⁱ
⟂R∥ = E∥ʳ / E_yⁱ
```
Both of these sets of equations are combined with the differential equation to be satisfied
by the wave fields at oblique incidence, giving two independent sets of relationships.

For **horizontal polarization**:
```math
(-A/C)*h₁(ζ) + (-B/C)*h₂(ζ) + E_yⁱ\\left(\\frac{⟂R⟂ + 1}{C}\\right) + E∥ⁱ \\left( \\frac{∥R⟂}{C} \\right) = 0
(-A/C)*(ch₁(ζ) + Kh₁′(ζ)) + (-B/C)*(ch₂(ζ) + Kh₂′(ζ)) + 2E_y = 0
```
For **vertical polarization**:
```math
(-Q/C)h₁(ζ) + (-G/C)h₂(ζ) + E∥ⁱ\\left(\\frac{∥R∥ + 1}{C}\\right) + E_yⁱ\\left(\\frac{⟂R∥}{C}\\right) = 0
(-Q/C)(ch₁(ζ) + (Kh₁′(ζ) + Lh₁(ζ))/n²) + (-G/C)(ch₂(ζ) + (Kh₂′(ζ) + Lh₂(ζ))/n²) + 2E∥ⁱ = 0
```
where we are to determine the coefficients `A`, `B`, `Q`, and `G`.
"""
function integratethroughfreespace(θ, k, z₀, zz, h, X₀)
    abs(zz-z₀) < 0.001 && return

    # TODO: Change z₀ and zz variable names

    # Initialize
    C = cosd(θ)
    C² = C^2
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2K)
    n₀² = 1 + α*(z₀ - h)
    nₜ² = 1 + α*(zz - h)

    # Complex "heights"
    ζ₀ = (k/α)^(2/3)*(C² + α*(z₀ - h))
    ζₜ = (k/α)^(2/3)*(C² + α*(zz - h))

    # (real(ζ₀) <= 0) & (real(ζₜ) <= 0) ? (bothnegative = true) : (bothnegative = false)

    # Computation of height-gain coefficients for two conditions on the upgoing wave from
    # height `z₀`; E∥ = 1, Ey = 0 and E∥ = 0, Ey = 1
    # mh1, mh2, mh1p, mh2p = modhankel(ζ₀, bothnegative)
    mh1 = modhankel1(ζ₀)
    mh2 = modhankel2(ζ₀)
    mh1p = modhankel1prime(ζ₀)
    mh2p = modhankel2prime(ζ₀)
    a₁ = @SMatrix [mh1            mh2;
                   C*mh1+K*mh1p   C*mh2+K*mh2p]
    Δ₁ = det(a₁)

    AC1 = X₀[2,1]*a₁[2,2]/Δ₁
    BC1 = -X₀[2,1]*a₁[2,1]/Δ₁
    AC2 = (X₀[2,2]*a₁[2,2] - 2a₁[1,2])/Δ₁
    BC2 = (2a₁[1,1] - X₀[2,2]*a₁[2,1])/Δ₁

    a₂ = @SMatrix [mh1                          mh2;
                   C*mh1+(K*mh1p+L*mh1)/n₀²     C*mh2+(K*mh2p+L*mh2)/n₀²]
    Δ₂ = det(a₂)

    QC1 = (X₀[1,1]*a₂[2,2] - 2a₂[1,2])/Δ₂
    GC1 = (2a₂[1,1] - X₀[1,1]*a₂[2,1])/Δ₂
    QC2 = X₀[1,2]*a₂[2,2]/Δ₂
    GC2 = -X₀[1,2]*a₂[2,1]/Δ₂

    # Computation of upgoing fields E∥ and Ey at height `zz` for the two conditions above.
    # mh1, mh2, mh1p, mh2p = modhankel(ζₜ, bothnegative)
    mh1 = modhankel1(ζₜ)
    mh2 = modhankel2(ζₜ)
    mh1p = modhankel1prime(ζₜ)
    mh2p = modhankel2prime(ζₜ)

    a21 = C*mh1 + K*mh1p
    a22 = C*mh2 + K*mh2p

    # Calculate parallel (p) and y fields
    # TODO: May be able to make some of this more matrix-y
    Eyt₁ = (AC1*a21 + BC1*a22)/2
    Eyt₂ = (AC2*a21 + BC2*a22)/2

    a21 = C*mh1 + (K*mh1p + L*mh1)/nₜ²
    a22 = C*mh2 + (K*mh2p + L*mh2)/nₜ²

    Ept₁ = (QC1*a21 + GC1*a22)/2
    Ept₂ = (QC2*a21 + GC2*a22)/2

    # Reflection matrix at the `zz` level
    W = Ept₁*Eyt₂ - Ept₂*Eyt₁
    V₁h = AC1*mh1 + BC1*mh2
    V₂h = AC2*mh1 + BC2*mh2
    V₁v = QC1*mh1 + GC1*mh2
    V₂v = QC2*mh1 + GC2*mh2

    Xₜ = @SMatrix [(V₁v*Eyt₂ - V₂v*Eyt₁)/W  (V₂v*Ept₁ - V₁v*Ept₂)/W;
                  (V₁h*Eyt₂ - V₂h*Eyt₁)/W   (V₂h*Ept₁ - V₁h*Ept₂)/W]

    return Xₜ
end

"""
Determine magnetoionic eigenvectors.

See NOSC TR 1143 section II. This is a combination of 1143 notation and LWPC algorithm.

ζₚ is only updated when in center of search rectangle and needs to be continuously returned
and repassed to function (or just a mutable struct). It should be initialized at the beginning
of program to [Complex(0.0), Complex(0.0)]. Same with Γₚ.

We use the form
```math
E_{ud} = \\begin{pmatrix}
            \\mathcal{H}_{y1} & W + ζ \\
            W + ζ & E_{y2}
        \\end{pmatrix} \\left( -2(W + ζ) \\right)^{1/2}
```
where the sign of the square root representing `ζ` is chosen so that the sum `W + ζ` has the
larger magnitude in any one search rectangle. Note that this `E` is singular for values of `θ`
for which `ζ = 0`, that is, for which ordinary waves cannot be distinguished from extraordinary
waves. Although not singular, it is not well conditioned for large collision frequency when
the propagation is east-west at the equator.
"""
function eigenmatrix(lenset, θ, ω, heightinterp, height, ζₚ, Γₚ, drcs::DirectionCosines)
    # Unpack
    dcl = drcs.dcl
    dcm = drcs.dcm
    dcn = drcs.dcn
    G = drcs.G

    # Initialize
    C = cosd(θ)
    S = sind(θ)
    C² = C^2
    S² = S^2
    dcl² = dcl^2
    dcm² = dcm^2
    dcn² = dcn^2

    α = 2(heightinterp - height)/earthradius
    p² = 1 + α
    q² = C² + α
    q = sqrt(q²)
    real(C) < 0 && (q = -q)

    e = @MArray zeros(ComplexF64, 2,2,2)
    inve = @MArray zeros(ComplexF64, 2,2,2)
    updown = +1  # +1 for upgoing waves, -1 for downgoing
    for ud = 1:2
        A = dcl*S + updown*dcn*q
        B = -dcn*S + updown*dcl*q

        W = (dcm²*p² - B^2)*G

        Hy₁ = 2(A + dcm*B*G)*p²
        Ey₂ = 2(A - dcm*B*G)

        ζ = sqrt(((B^2 + dcm²*p²)*G)^2 - 4A^2*p²)

        # `lenset` is set to `true` in `introe()` for θ at the center of the search rectangle
        # and set to `false` for all other rectangle locations. The sign of `ζ` is chosen so
        # that `W + ζ` has the larger magnitude. The sign of `ζ` elsewhere is chosen so that
        # the value of `ζ` is closer to that of `ζ` at the center.
        lenset && (ζₚ[ud] = W)
        abs2(ζ-ζₚ[ud]) > abs2(ζ+ζₚ[ud]) && (ζ = -ζ)
        lenset && (ζₚ[ud] = ζ)

        # THe sign of Γ is chosen so that the value of Γ is closer to that of Γ at the center
        # of the rectangle
        Γ = sqrt(-2*ζ*(W + ζ))

        lenset && (Γₚ[ud] = Γ)
        abs2(Γ-Γₚ[ud]) > abs2(Γ+Γₚ[ud]) && (Γ = -Γ)

        # Eigenvector matrix
        e[1,1,ud] = Hy₁
        e[1,2,ud] = W+ζ
        e[2,1,ud] = W+ζ
        e[2,2,ud] = Ey₂
        e[:,:,ud] ./= Γ

        inve[:,:,ud] = inv(e[:,:,ud])  # TODO: Is this actually needed or can we eventually do a matrix division?

        # Change the sign for downgoing waves
        updown = -1
    end

    return ζₚ, Γₚ, e, inve
end

"""
Interpolates values of the elements of the magnetoionic reflection matrix.

See section V of NOSC TR 1143.
"""
function initoe()

end

end
