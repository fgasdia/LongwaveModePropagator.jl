using SpecialFunctions
using LinearAlgebra
using StaticArrays
using PolynomialRoots
using DifferentialEquations

# """
# LWPC's MF_WVGD is a medium length subroutine that obtains approximate eigenangles (_modes_)
# in a single path segment. It calls several secondary routines to accomplish this task.
# """
# function waveguide(inputs::Inputs)
# # Unpack
# azim, dip = inputs.azim, inputs.dip
#
# # Initialize
# dcl = cosd(dip)*cosd(azim)
# dcm = cosd(dip)*sind(azim)
# dcn = -sind(dip)
#
# # Initialize search area
# Zb, Ze = boundaries(freq)
# Δθmesh = sqrt(3.75e3/freq)  # (deg) search grid size
# tolerance = 0.1  # (deg)
# end

"""
Search boundary for zeros in complex plane.

TODO: This seems a little arbitrary...

See also: `lwp_input.for`
"""
function boundaries(freq)
    if freq < 20e3
        Zb = complex(60., 0.)
        Ze = complex(90., -9.)
        return Zb, Ze
    else
        Zb = complex(40., 0.)
        Ze = complex(90., -6.)
        return Zb, Ze
    end
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
function tmatrix(θ, ω, referenceheight, height, spec::Constituent, B, dcl, dcm, dcn)
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

    # In LWPC and Sheddy 1968 Fortran Program `earthcurvature` is not multiplied by capd
    # (the above line), even though it _is_ multiplied in MS 1976.
    # This seems to be supported by Pappert 1968
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

See LWPC: `mf_initr.for`
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
        angq = angle(qval)*180/π
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

    X = MMatrix{2, 2}(2(T[1]*(C + q[2]) - T[2]*(C +q[1]))/Δ,
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
    (-im*k/2)*(A + B*X + X*C + X*D*X)
end

"""
"""
function xderivative(X, params, height)
    T = tmatrix(params.θ, params.ω, params.referenceheight, height,
                params.species, params.B, params.dcl, params.dcm, params.dcn)[2]
    A, B, C, D = smatrix(params.θ, T)
    dXdh(X, A, B, C, D, params.k)
end

"""
Full wave solution for `X` through ionosphere.

Integrates ``dX/dh`` (`dXdh()`) from `topheight` to `bottomheight`.
"""
function integratethroughionosphere(θ, ω, k, fromheight, toheight, referenceheight,
                                    species, B, dcl, dcm, dcn)
    # Initial condition for integration
    # `fromheight` should be topheight
    M = tmatrix(θ, ω, referenceheight, fromheight, species, B, dcl, dcm, dcn)[1]
    initialX = sharplyboundedreflectionmatrix(θ, M)[2]

    params = (θ=θ, referenceheight=referenceheight, ω=ω, k=k, species=species,
              B=B, dcl=dcl, dcm=dcm, dcn=dcn)

    heightbounds = (fromheight, toheight)
    prob = ODEProblem(xderivative, initialX, heightbounds, params)
    # sol = solve(prob, BS3(), reltol=1e-3, dtmax=1., save_everystep=false, save_start=false)  # Faster, less accurate
    sol = solve(prob, Tsit5(), reltol=1e-4, save_everystep=false)  # Decent accuracy for the speed

    # TODO: Catch if `sol` is an error, then `nointegration = true`, else `nointegration = false`

    return sol
end
# function integratethroughionosphere(mode::Mode, inputs::Inputs, referenceheight,
#                                     species::Constituent, B, drcs::DirectionCosines)
#     # Unpack
#     θ, ω, k = mode.θ, mode.ω, mode.wavenumber
#     bottomheight, topheight = inputs.bottomheight, inputs.topheight
#     dcl, dcm, dcn = drcs.dcl, drcs.dcm, drcs.dcn
#
#     integratethroughionosphere(θ, ω, k, topheight, bottomheight, species, B, dcl, dcm, dcn)
# end

"""
These functions calculate the modified Hankel functions of order 1/3 and their derivatives.

The functions h₁(ζ) and h₂(ζ) satisfy the Stokes differential equation (_Airy_ function)
```math
\\frac{d²w}{dζ²} + ζw = 0
```
"""
function modhankel(z)
    Ai = SpecialFunctions.airyai(-z)
    Bi = SpecialFunctions.airybi(-z)
    Ai′ = SpecialFunctions.airyaiprime(-z)
    Bi′ = SpecialFunctions.airybiprime(-z)

    mh1 = 12^(1/6)*exp(-im*π/6)*(Ai-im*Bi)
    mh2 = 12^(1/6)*exp(im*π/6)*(Ai+im*Bi)
    mh1p = Complex(-1)^(5/6)*12^(1/6)*(Ai′-im*Bi′)
    mh2p = Complex(-1)^(-5/6)*12^(1/6)*(Ai′+im*Bi′)

    return mh1, mh2, mh1p, mh2p
end
modhankel1(z) = 12^(1/6)*exp(-im*π/6)*(SpecialFunctions.airyai(-z)-im*SpecialFunctions.airybi(-z))
modhankel2(z) = 12^(1/6)*exp(im*π/6)*(SpecialFunctions.airyai(-z)+im*SpecialFunctions.airybi(-z))
modhankel1prime(z) = Complex(-1)^(5/6)*12^(1/6)*(SpecialFunctions.airyaiprime(-z)-im*SpecialFunctions.airybiprime(-z))
modhankel2prime(z) = Complex(-1)^(-5/6)*12^(1/6)*(SpecialFunctions.airyaiprime(-z)+im*SpecialFunctions.airybiprime(-z))

"""
Performs an integration of the differential equations for the ionosphere reflection matrix
through a free space region over a curved earth.

The integration may be performed in either a positive or negative height direction, but LWPC
always integrates upward. Therefore, argument `X` is actually `X₀`, the modified reflection
coefficient matrix at `bottomheight`. The integration variables are the matrix ``X = (R+1)/C``
where ``R`` is the reflection matrix described by Budden. This solution is based on Budden,
Radio Waves in the Ionosphere, in particular the material on pg 118, 327-329, 336-338, and
343-345. MS 1976, especially Appendix 2, presents a detailed explanation.

Computation at θ = 90° is not excluded. Also, the equations are formulated so that a smooth
transition is made from the "curved earth" form to the "flat earth" form for appropriate
values of θ.

TODO: Is this step even needed? By integrating up to an effective reflection height, we
guarantee no poles in function, but can we just apply more pure reflection matrix in modal
function since we have a root finder that can handle poles?

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
function integratethroughfreespace!(X, θ, k, fromheight, toheight, referenceheight)
    abs(toheight - fromheight) < 0.001 && return X

    # Initialize
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2k)

    C = cosd(θ)
    C² = C^2

    # Complex "heights"
    ζ₀ = (k/α)^(2/3)*(C² + α*(fromheight - referenceheight))
    n₀² = 1 + α*(fromheight - referenceheight)

    # Computation of height-gain coefficients for two conditions on the upgoing wave from
    # `fromheight`; E∥ = 1, Ey = 0 and E∥ = 0, Ey = 1
    mh1, mh2, mh1p, mh2p = modhankel(ζ₀)

    a₁ = @SMatrix [mh1            mh2;
                   C*mh1+K*mh1p   C*mh2+K*mh2p]
    Δ₁ = det(a₁)

    AC1 = X[2,1]*a₁[2,2]/Δ₁
    BC1 = -X[2,1]*a₁[2,1]/Δ₁
    AC2 = (X[2,2]*a₁[2,2] - 2*a₁[1,2])/Δ₁
    BC2 = (2*a₁[1,1] - X[2,2]*a₁[2,1])/Δ₁

    a₂ = @SMatrix [mh1                          mh2;
                   C*mh1+(K*mh1p+L*mh1)/n₀²     C*mh2+(K*mh2p+L*mh2)/n₀²]
    Δ₂ = det(a₂)

    QC1 = (X[1,1]*a₂[2,2] - 2*a₂[1,2])/Δ₂
    GC1 = (2*a₂[1,1] - X[1,1]*a₂[2,1])/Δ₂
    QC2 = X[1,2]*a₂[2,2]/Δ₂
    GC2 = -X[1,2]*a₂[2,1]/Δ₂

    # Computation of upgoing fields E∥ and Ey at `toheight` for the two conditions above.
    ζₜ = (k/α)^(2/3)*(C² + α*(toheight - referenceheight))
    nₜ² = 1 + α*(toheight - referenceheight)

    mh1, mh2, mh1p, mh2p = modhankel(ζₜ)

    a21 = C*mh1 + K*mh1p
    a22 = C*mh2 + K*mh2p

    # Calculate parallel (p) and y fields
    Eyt₁ = (AC1*a21 + BC1*a22)/2
    Eyt₂ = (AC2*a21 + BC2*a22)/2

    a21 = C*mh1 + (K*mh1p + L*mh1)/nₜ²
    a22 = C*mh2 + (K*mh2p + L*mh2)/nₜ²

    Ept₁ = (QC1*a21 + GC1*a22)/2
    Ept₂ = (QC2*a21 + GC2*a22)/2

    # Reflection matrix at the `toheight` level
    W = Ept₁*Eyt₂ - Ept₂*Eyt₁
    V₁h = AC1*mh1 + BC1*mh2
    V₂h = AC2*mh1 + BC2*mh2
    V₁v = QC1*mh1 + GC1*mh2
    V₂v = QC2*mh1 + GC2*mh2

    # Mutate to `X` at the `toheight` level
    X[1,1] = V₁v*Eyt₂ - V₂v*Eyt₁
    X[1,2] = V₂v*Ept₁ - V₁v*Ept₂
    X[2,1] = V₁h*Eyt₂ - V₂h*Eyt₁
    X[2,2] = V₂h*Ept₁ - V₁h*Ept₂
    X ./= W  # XXX: Is there really just 1 W?

    return X
end
# function integratethroughfreespace!()
#     # Unpack
#     θ, k = mode.θ, mode.wavenumber
#
#     integratethroughfreespace!()
# end

r2x(R,C) = (R+I)/C
x2r(X,C) = C*X - I

"""
Calculate `n` and `d`, the modified ground reflection matrix, at `referenceheight`.

This is based on the Fresnel reflection coefficients for the ground/free-space interface.
The modified matrix terms are given by:
```nd⁻¹ = (Rg⁻¹ + I)/C```
or
```
∥n∥/∥d∥ = \\frac{∥R̄∥ + 1}{∥R̄∥C}
```
and
```
⟂n⟂/⟂d⟂ = \\frac{⟂R̄⟂ + 1}{⟂R̄⟂C}
```
The ground reflection matrix is therefore ``R̄11 = D11/(cosθ*N11 - D11)``.

Morfitt Shellman 1976, pg 25, and Appendix A and C includes information on the derivation.
Also look at Pappert et al 1966, "A Numerical Investigation of Classical Approximations used
in VLF Propagation". Much of this code does not follow the math in MS 1976, but is based
upon `mf_bars.for`.

`toheight` should be `reflectionheight` in normal usage

See also: `mf_rbars.for`
"""
function groundreflection(θ, ω, k, σ, ϵᵣ, toheight, referenceheight)
    # Initialize
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2k)

    C = cosd(θ)
    C² = C^2

    n11, n22, d11, d22 = fresnelgroundreflection(θ, ω, σ, ϵᵣ)

    ζ₀ = (k/α)^(2/3)*(C² - α*referenceheight)
    n₀² = 1 - α*referenceheight
    a₀ = C + L/n₀²
    h1₀, h2₀, h1₀′, h2₀′ = modhankel(ζ₀)

    ζₜ = (k/α)^(2/3)*(C² + α*(toheight - referenceheight))
    nₜ² = 1 + α*(toheight - referenceheight)
    aₜ = C + L/nₜ²
    h1ₜ, h2ₜ, h1ₜ′, h2ₜ′ = modhankel(ζₜ)

    k1₀ = K/n₀²
    k1ₜ = K/nₜ²

    f = h1₀*h2ₜ - h2₀*h1ₜ
    f0 = h1₀′*h2ₜ - h2₀′*h1ₜ
    ft = h1₀*h2ₜ′ - h2₀*h1ₜ′
    f0t = h1₀′*h2ₜ′ - h2₀′*h1ₜ′

    num11 = -2*n11*(a₀*f + k1₀*f0) + 4*d11*f
    num22 = -2*n22*(C*f + K*f0) + 4*d22*f
    den11 = -n11*(a₀*aₜ*f + k1ₜ*a₀*ft + k1₀*aₜ*f0 + k1₀*k1ₜ*f0t) + 2*d11*(aₜ*f + k1ₜ*ft)
    den22 = -n22*(C²*f + K*C*ft + K*C*f0 + K^2*f0t) + 2*d22*(C*f + K*ft)

    return num11, num22, den11, den22
end

"""
Fresnel reflection matrix at the ground in terms of `n` and `d`.

See Morfitt Shellman 1976, Appendix C for the derivation.
"""
function fresnelgroundreflection(θ, ω, σ, ϵᵣ)
    # Initialize
    C = cosd(θ)
    S = sind(θ)
    C² = C^2
    S² = S^2

    # At the "from" level, ``z = 0``
    ng² = complex(ϵᵣ, -σ/(ω*ϵ₀))
    W = sqrt(ng² - S²)

    n11 = 1
    n22 = 1/W
    d11 = 0.5*(C - W/ng²)
    d22 = 0.5*(C/W - 1)

    return n11, n22, d11, d22
end

"""
Note: the values `N11`, `N22`, `D11`, and `D22` may not match those returned from LWPC's
`mf_rbars.for` subroutine, but the ratios ``N11/D11`` and ``N22/D22`` are both equal.
"""
function groundreflectionms76(θ, ω, k, σ, ϵᵣ, toheight, referenceheight)
    # Initialize
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2k)

    C = cosd(θ)
    S = sind(θ)
    C² = C^2
    S² = S^2

    # At the "from" level, ``z = 0``
    ng² = complex(ϵᵣ, -σ/(ω*ϵ₀))
    W = sqrt(ng² - S²)

    n11 = 1
    n22 = 1/W
    d11 = 0.5*(C - W/ng²)
    d22 = 0.5*(C/W - 1)

    ζ₀ = (k/α)^(2/3)*(C² - α*referenceheight)
    n₀² = 1 - α*referenceheight

    mh1, mh2, mh1p, mh2p = modhankel(ζ₀)
    e = -im*(2/3)*ζ₀^(3/2) + im*π*5/12
    mh1 /= exp(-e)
    mh2 /= exp(e)
    mh1p /= exp(-e)
    mh2p /= exp(e)

    # For horizontal polarization
    a₁ = @SMatrix [mh1              mh2;
                   C*mh1+K*mh1p     C*mh2+K*mh2p]
    Δ₁ = det(a₁)

    AC2 = (n22*a₁[2,2] - 2*d22*a₁[1,2])/(Δ₁*d22)
    BC2 = (2*d22*a₁[1,1] - n22*a₁[2,1])/(Δ₁*d22)

    # For vertical polarization
    a₂ = @SMatrix [mh1                          mh2;
                   C*mh1+(K*mh1p+L*mh1)/n₀²     C*mh2+(K*mh2p+L*mh2)/n₀²]
    Δ₂ = det(a₂)

    QC1 = (n11*a₂[2,2] - 2*d11*a₂[1,2])/(Δ₂*d11)
    GC1 = (2*d11*a₂[1,1] - n11*a₂[2,1])/(Δ₂*d11)

    # At the "to" level, ``z = toheight``
    ζₜ = (k/α)^(2/3)*(C² + α*(toheight - referenceheight))
    nₜ² = 1 + α*(toheight - referenceheight)

    mh1, mh2, mh1p, mh2p = modhankel(ζₜ)
    e = -im*(2/3)*ζₜ^(3/2) + im*π*5/12
    mh1 /= exp(-e)
    mh2 /= exp(e)
    mh1p /= exp(-e)
    mh2p /= exp(e)

    # Horizontal polarization
    a21 = C*mh1 + K*mh1p
    a22 = C*mh2 + K*mh2p

    Eyt₂ = (AC2*a21 + BC2*a22)/2

    # Vertical polarization
    a21 = C*mh1 + (K*mh1p + L*mh1)/nₜ²
    a22 = C*mh2 + (K*mh2p + L*mh2)/nₜ²

    Ept₁ = (QC1*a21 + GC1*a22)/2

    N11 = QC1*mh1 + GC1*mh2
    N22 = AC2*mh1 + BC2*mh2
    D11 = Ept₁
    D22 = Eyt₂

    return N11, N22, D11, D22
end

"""
Modified modal function "F₃(θ)" that has no poles and same zeros as physical modal equation except
no zeros at ``θ = 90°``.

A full description of this function and its derivation is given in Morfitt Shellman 1976
(although note that F₃(θ) given at the beginning of the report is incorrect, see Appendix 1)
```
F₃(θ) = (∥n∥ - ∥X∥∥d∥)(⟂n⟂ - ⟂X⟂⟂d⟂) - ∥X⟂⟂X∥∥d∥⟂d⟂
```
where ``X = (R+I)/c``, ``c=cos(θ)``, and ``nd⁻¹ = (Rg⁻¹ + I)/c``.

TODO: Follow ELF approach, except for VLF/LF. Integrate through ionosphere down to the ground.
No further integration back up to a reference height is done. Zeros are found at z=0.

See also: `mf_fctval.for`, `groundreflection`(@ref)
"""
function modifiedmodalfunction(θ, ω, k, σ, ϵᵣ,
                               bottomheight, topheight, reflectionheight, referenceheight,
                               species,
                               B, dcl, dcm, dcn)

    # Calculate `X` at `bottomheight`
    Xsol = integratethroughionosphere(θ, ω, k, topheight, bottomheight, referenceheight,
                                      species, B, dcl, dcm, dcn)
    X = Xsol[end]

    # Refer `X` from `bottomheight` to `referenceheight`
    integratethroughfreespace!(X, θ, k, bottomheight, reflectionheight, referenceheight)

    # Calculate ground reflection matrix as `N` and `D` referred to `referenceheight`
    N11, N22, D11, D22 = groundreflection(θ, ω, k, σ, ϵᵣ, reflectionheight, referenceheight)

    f = (N11 - X[1,1]*D11)*(N22 - X[2,2]*D22) - X[2,1]*X[1,2]*D11*D22
end

"""
Calculate modal excitation factors.

TODO: Support horizontal end and broadside exciters (antenna).

See Morfitt Shellman 1976 pg 29.
"""
function excitationfactors()
    # Initialize
    S = sind(θ)

    B₁ = S^(5/2)/dFdθ
    B₂ = -B₁/S

    Ez = B₁*(1+Rg11)^2*(1 - Rg22*R22)/(Rg11*D11)
    Ey = -B₁/S * R21*(1 + Rg11)*(1 + Rg22)/D12
    Ex = B₁/S * (1 + Rg11)^2*(1 - R22^2)/(Rg11*D11)

    return Ez, Ey, Ex
end

"""
Calculate modal height gain terms.
"""
function heightgain()
    ζ = (k/α)^(2/3)*(C² + α*(Z - H))

    mh1 = modhankel1(ζ)
    mh2 = modhankel2(ζ)
    fvertical(z) = (QC1*mh1 + GC*mh2)*exp((z - d)/earthradius)
    fhorizontal(z) = (AC2*mh1 + BC2*mh2)
    g(z) = -im*exp((z - d)/a)*(cbrt(2/(a*k))*(QC1*mh1p + GC1*mh2p) + (2/(a*k))*(QC1*mh1 + GC1*mh2))
end
