using SpecialFunctions
using LinearAlgebra
using StaticArrays
using PolynomialRoots
using DifferentialEquations
using Markdown

using ModifiedHankelFunctionsOfOrderOneThird

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

struct GroundCoefficients{S,T}
    C::T
    W::T
    ng²::T
    K::T
    L::T
    cbrtkoverαsq::S
    ζ₀::T
    ζₜ::T
    n₀²::S
    nₜ²::S
    B1₀::T
    B2₀::T
    B3₀::T
    B4₀::T
    A1::T
    A2::T
    A3::T
    A4::T
    B1ₜ::T
    B2ₜ::T
    B3ₜ::T
    B4ₜ::T
    h1₀::T
    h2₀::T
    h1₀′::T
    h2₀′::T
    h1ₜ::T
    h2ₜ::T
    h1ₜ′::T
    h2ₜ′::T
end

"""
Search boundary for zeros in complex plane.

TODO: This seems a little arbitrary...

See also: `lwp_input.for`
"""
function boundaries(freq)
    if freq < 20e3
        Zb = complex(60, 0)
        Ze = complex(90, -9)
        return Zb, Ze
    else
        Zb = complex(40, 0)
        Ze = complex(90, -6)
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
"""
function mmatrix(ω, referenceheight, height, spec::Constituent, B, dcl, dcm, dcn)
    # Initialize
    earthcurvature = 2/earthradius*(height - referenceheight)

    # TODO: Add M up for each species
    # Constitutive relations (see Budden 1955, pg. 517)
    e, m, N, ν = spec.charge, spec.mass, spec.numberdensity, spec.collisionfrequency
    X = N(height)*e^2/(ϵ₀*m*ω^2)
    Y = e*B/(m*ω)
    Z = ν(height)/ω
    U = 1 - im*Z

    # Construct susceptibility matrix (see Budden 1955, eq. 3)
    U² = U^2
    Y² = Y^2

    M = @MMatrix [U²-dcl^2*Y²               -im*dcn*Y*U-dcl*dcm*Y²    im*dcm*Y*U-dcl*dcn*Y²;
                  im*dcn*Y*U-dcl*dcm*Y²     U²-dcm^2*Y²               -im*dcl*Y*U-dcm*dcn*Y²;
                  -im*dcm*Y*U-dcl*dcn*Y²    im*dcl*Y*U-dcm*dcn*Y²     U²-dcn^2*Y²]

    M .*= -X/(U*(U² - Y²))

    # In LWPC and Sheddy 1968 Fortran Program `earthcurvature` is not multiplied by capd
    # (the above line), even though it _is_ multiplied in MS 1976.
    # This seems to be supported by Pappert 1968
    M[1,1] += earthcurvature
    M[2,2] += earthcurvature
    M[3,3] += earthcurvature

    return M
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

Matrix ``T`` fulfills the differential equation ``e′ = -iTe`` where ``e`` is the column
matrix ``e = [Eₓ -Ey Hₓ Hy]``.

Note:
- `T44` in MS 76 should be negative
- MS 76' says ``C = -S^{11} + S^{12}`` and then writes that in terms of ``T``, but `C11` is
missing a minus sign, it should be ``-T^{11} - CT^{41}``.
"""
function smatrix(θ, M)
    c = cosd(θ)
    s = sind(θ)
    c² = c^2

    oneplusM33 = 1 + M[3,3]

    T11 = -s*M[3,1]/oneplusM33
    T12 = s*M[3,2]/oneplusM33
    # T13 = 0
    T14 = (c² + M[3,3])/oneplusM33
    # T21 = 0
    # T22 = 0
    # T23 = 1
    # T24 = 0
    T31 = M[2,3]*M[3,1]/oneplusM33 - M[2,1]
    T32 = c² + M[2,2] - M[2,3]*M[3,2]/oneplusM33
    # T33 = 0
    T34 = s*M[2,3]/oneplusM33
    T41 = 1 + M[1,1] - M[1,3]*M[3,1]/oneplusM33
    T42 = M[3,2]*M[1,3]/oneplusM33 - M[1,2]
    # T43 = 0
    T44 = -s*M[1,3]/oneplusM33

    A = @SMatrix [4T41   0;
                  0      4]
    B = @SMatrix [2(T44-c*T41)    -2T42;
                  0               -2c]
    C = @SMatrix [-2(T11+c*T41)   0;
                  2T31            -2c]
    D = @SMatrix [c*(T11-T44)-T14+c²*T41   T12+c*T42;
                  -c*T31+T34               c²-T32]

    return A, B, C, D
end

function dsmatrixdC(θ, M)
    c = cosd(θ)
    s = sind(θ)

    ds = -c/s  # ds/dc
    dcs = s - c^2/s

    oneplusM33 = 1 + M[3,3]

    T31 = M[2,3]*M[3,1]/oneplusM33 - M[2,1]
    T41 = 1 + M[1,1] - M[1,3]*M[3,1]/oneplusM33

    dA = @SMatrix zeros(eltype(M), 2, 2)
    dB = @SMatrix [-2(T41+M[1,3]*ds/oneplusM33)   0;
                   0                              -2]
    dC = @SMatrix [-2(T41-M[3,1]*ds/oneplusM33)   0;
                   0                              -2]
    dD = @SMatrix [2*c*(T41-1+M[3,3]/oneplusM33)+dcs*(-M[3,1]+M[1,3])/oneplusM33    ds*M[3,2]/oneplusM33;
                   -T31+ds*M[2,3]/oneplusM33                                        0]

    return dA, dB, dC, dD
end


"""
Calculate angle from 315°.

Used to sort in [`sharplyboundedreflectionmatrix`](@ref).
"""
function anglefrom315(qval)
    angq = rad2deg(angle(qval))
    angq < 0 && (angq += 360)
    angq < 135 && (angq += 360)
    abs(angq - 315)
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
    C² = C^2
    P = MVector{2,ComplexF64}(undef)
    T = MVector{2,ComplexF64}(undef)

    # Booker quartic coefficients
    b4 = 1 + M[3,3]
    b3 = S*(M[1,3] + M[3,1])
    b2 = -(C² + M[3,3])*(1 + M[1,1]) + M[1,3]*M[3,1] -
         (1 + M[3,3])*(C² + M[2,2]) + M[2,3]*M[3,2]
    b1 = S*(M[1,2]*M[2,3] + M[2,1]*M[3,2] -
         (C² + M[2,2])*(M[1,3] + M[3,1]))
    b0 = (1 + M[1,1])*(C² + M[2,2])*(C² + M[3,3]) +
         M[1,2]*M[2,3]*M[3,1] + M[1,3]*M[2,1]*M[3,2] -
         M[1,3]*(C² + M[2,2])*M[3,1] -
         (1 + M[1,1])*M[2,3]*M[3,2] -
         M[1,2]*M[2,1]*(C² + M[3,3])

    q = roots([b0, b1, b2, b3, b4])
    sort!(q, by=anglefrom315)

    # Dispersion matrix elements for R, X
    D12 = M[1,2]
    D32 = M[3,2]
    D33 = C² + M[3,3]
    @inbounds for j = 1:2
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
    -im*k/2*(A + B*X + X*C + X*D*X)
end

"""
"""
function xderivative(X, params, height)
    M = mmatrix(params.ω, params.referenceheight, height,
                params.species, params.B, params.dcl, params.dcm, params.dcn)
    A, B, C, D = smatrix(params.θ, M)
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
    M = mmatrix(ω, referenceheight, fromheight, species, B, dcl, dcm, dcn)
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
    h1, h2, h1′, h2′ = modifiedhankel(ζ₀)

    a₁ = @SMatrix [h1            h2;
                   C*h1+K*h1′    C*h2+K*h2′]
    Δ₁ = det(a₁)

    AC1 = X[2,1]*a₁[2,2]/Δ₁
    BC1 = -X[2,1]*a₁[2,1]/Δ₁
    AC2 = (X[2,2]*a₁[2,2] - 2*a₁[1,2])/Δ₁
    BC2 = (2*a₁[1,1] - X[2,2]*a₁[2,1])/Δ₁

    a₂ = @SMatrix [h1                          h2;
                   C*h1+(K*h1′+L*h1)/n₀²       C*h2+(K*h2′+L*h2)/n₀²]
    Δ₂ = det(a₂)

    QC1 = (X[1,1]*a₂[2,2] - 2*a₂[1,2])/Δ₂
    GC1 = (2*a₂[1,1] - X[1,1]*a₂[2,1])/Δ₂
    QC2 = X[1,2]*a₂[2,2]/Δ₂
    GC2 = -X[1,2]*a₂[2,1]/Δ₂

    # Computation of upgoing fields E∥ and Ey at `toheight` for the two conditions above.
    ζₜ = (k/α)^(2/3)*(C² + α*(toheight - referenceheight))
    nₜ² = 1 + α*(toheight - referenceheight)

    h1, h2, h1′, h2′ = modifiedhankel(ζₜ)

    a21 = C*h1 + K*h1′
    a22 = C*h2 + K*h2′

    # Calculate parallel (p) and y fields
    Eyt₁ = (AC1*a21 + BC1*a22)/2
    Eyt₂ = (AC2*a21 + BC2*a22)/2

    a21 = C*h1 + (K*h1′ + L*h1)/nₜ²
    a22 = C*h2 + (K*h2′ + L*h2)/nₜ²

    Ept₁ = (QC1*a21 + GC1*a22)/2
    Ept₂ = (QC2*a21 + GC2*a22)/2

    # Reflection matrix at the `toheight` level
    W = Ept₁*Eyt₂ - Ept₂*Eyt₁
    V₁h = AC1*h1 + BC1*h2
    V₂h = AC2*h1 + BC2*h2
    V₁v = QC1*h1 + GC1*h2
    V₂v = QC2*h1 + GC2*h2

    # Mutate to `X` at the `toheight` level
    X[1,1] = V₁v*Eyt₂ - V₂v*Eyt₁
    X[1,2] = V₂v*Ept₁ - V₁v*Ept₂
    X[2,1] = V₁h*Eyt₂ - V₂h*Eyt₁
    X[2,2] = V₂h*Ept₁ - V₁h*Ept₂
    X ./= W  # XXX: Is there really just 1 W?

    return X
end

r2x(R,C) = (R+I)/C
x2r(X,C) = C*X - I

"""
    _groundcoeffs(θ, ω, k, σ, ϵᵣ, toheight, referenceheight)

Calculate values needed for calculating `n` and `d`, the modified ground reflection matrix
terms, and their derivatives.

To convert between the coefficients used in this code and used in the MS76 derivation:

| B₁ | a₂(2,1)      | Ch₁₀ + (Kh₁₀′ + Lh₁₀)/n₀²  |
| B₂ | a₂(2,2)      | Ch₂₀ + (Kh₂₀′ + Lh₂₀)/n₀²  |
| B₃ | a₁(2,1)      | Ch₁₀ + Kh₁₀′               |
| B₄ | a₁(2,2)      | Ch₂₀ + Kh₂₀′               |
| A₁ | -GC1 (Δ₂d11) | n11 a₂(2,1) - 2d11 a₂(1,1) |
| A₂ | QC2 (Δ₂d11)  | n11 a₂(2,2) - 2d11 a₂(1,2) |
| A₃ | -BC2 (Δ₁d22) | n22 a₁(2,1) - 2d22 a₁(1,1) |
| A₄ | AC2 (Δ₁d22)  | n22 a₁(2,2) - 2d22 a₁(1,2) |
| B₁ | vert a21     | Ch₁ₜ + (Kh₁ₜ′ + Lh₁ₜ)/nₜ²     |
| B₂ | vert a22     | Ch₂ₜ + (Kh₂ₜ′ + Lh₂ₜ)/nₜ²     |
| B₃ | hor a21      | Ch₁ₜ + Kh₁ₜ′                 |
| B₄ | hor a22      | Ch₂ₜ + Kh₂ₜ′                 |

`toheight` should be `reflectionheight` in normal usage

References:
    - Morfitt Shellman 1976, pg. 25, Appendix A, Appendix C, pg. 195
    - Pappert et al 1966, "A Numerical Investigation of Classical Appoximations..."

See also: `mf_rbars.for`, [`ndground`](@ref)
"""
function _groundcoeffs(θ, ω, k, σ, ϵᵣ, toheight, referenceheight)
    # Initialize
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2k)
    cbrtkoverαsq = cbrt(k/α)^2

    # Initialize
    C = cosd(θ)
    S = sind(θ)
    C² = C^2

    # At the "from" level, ``z = 0`` (ground)
    ng² = complex(ϵᵣ, -σ/(ω*ϵ₀))
    W = sqrt(ng² - S^2)

    n11, n22, d11, d22 = fresnelnd(C, W, ng²)

    # At the ground
    ζ₀ = cbrtkoverαsq*(C² - α*referenceheight)
    n₀² = 1 - α*referenceheight
    h1₀, h2₀, h1₀′, h2₀′ = modifiedhankel(ζ₀)

    # At the "to" level, ``z = toheight``
    ζₜ = cbrtkoverαsq*(C² + α*(toheight - referenceheight))
    nₜ² = 1 + α*(toheight - referenceheight)
    h1ₜ, h2ₜ, h1ₜ′, h2ₜ′ = modifiedhankel(ζₜ)

    B1₀ = C*h1₀+(K*h1₀′+L*h1₀)/n₀²
    B2₀ = C*h2₀+(K*h2₀′+L*h2₀)/n₀²
    B3₀ = C*h1₀+K*h1₀′
    B4₀ = C*h2₀+K*h2₀′

    A1 = n11*B1₀ - 2*d11*h1₀
    A2 = n11*B2₀ - 2*d11*h2₀
    A3 = n22*B3₀ - 2*d22*h1₀
    A4 = n22*B4₀ - 2*d22*h2₀

    B1ₜ = C*h1ₜ + (K*h1ₜ′ + L*h1ₜ)/nₜ²
    B2ₜ = C*h2ₜ + (K*h2ₜ′ + L*h2ₜ)/nₜ²
    B3ₜ = C*h1ₜ + K*h1ₜ′
    B4ₜ = C*h2ₜ + K*h2ₜ′

    return GroundCoefficients(C, W, ng², K, L, cbrtkoverαsq, ζ₀, ζₜ, n₀², nₜ²,
                              B1₀, B2₀, B3₀, B4₀, A1, A2, A3, A4, B1ₜ, B2ₜ, B3ₜ, B4ₜ,
                              h1₀, h2₀, h1₀′, h2₀′, h1ₜ, h2ₜ, h1ₜ′, h2ₜ′)
end

doc"""
    ndground(A1, A2, A3, A4, B1, B2, B3, B4, h1ₜ, h2ₜ)

Calculate `n` and `d`, the modified ground reflection matrix terms, at `referenceheight`,
given the common ground coefficients.

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
The ground reflection matrix is therefore ``R̄11 = D11/(cos(θ)*N11 - D11)``.

From MS1976 we can derive
```math
_\parallel X_g_\parallel = \frac{\texttt{QC}_1 h_1 + \texttt{GC}_1 h_2}{\texttt{E}_\parallel_t_1}
_\perp X_g_\perp = \frac{\texttt{AC}_2 h_1 + \texttt{BC}_2 h_2}{E_y_t_2}
```

The values for `N11`, `N22`, `D11`, and `D22` found here match LWPC, but do not match the
results of the math derived in MS76. The ratios ``N11/D11`` and ``N22/D22`` remain the same.
    - `N11` is 2Δ₂d11 times derived N11
    - `N22` is 2Δ₁d22 times derived N22
    - `D11` is 2Δ₂d11 times derived D11
    - `D22` is 2Δ₁d22 times derived D22

References:
    - Morfitt Shellman 1976, pg. 25, Appendix A, Appendix C, pg. 195
    - Pappert et al 1966, "A Numerical Investigation of Classical Appoximations..."

See also: `mf_rbars.for`, [`_groundcoeffs`](@ref)
"""
function ndground(coeffs::GroundCoefficients)
    A1, A2, A3, A4 = coeffs.A1, coeffs.A2, coeffs.A3, coeffs.A4
    B1ₜ, B2ₜ, B3ₜ, B4ₜ = coeffs.B1ₜ, coeffs.B2ₜ, coeffs.B3ₜ, coeffs.B4ₜ
    h1ₜ, h2ₜ = coeffs.h1ₜ, coeffs.h2ₜ

    N11 = 2(A2*h1ₜ - A1*h2ₜ)
    N22 = 2(A4*h1ₜ - A3*h2ₜ)
    D11 = A2*B1ₜ - A1*B2ₜ
    D22 = A4*B3ₜ - A3*B4ₜ

    return N11, N22, D11, D22
end

"""
Return the derivative of the ground reflection coefficient terms `n` and `d` wrt θ at the
reflection height.

See also: [`ndground`](@ref)
"""
function dndgrounddC(coeffs::GroundCoefficients)
    # Unpack
    C, W, ng² = coeffs.C, coeffs.W, coeffs.ng²
    K, L, cbrtkoverαsq = coeffs.K, coeffs.L, coeffs.cbrtkoverαsq
    ζ₀, ζₜ = coeffs.ζ₀, coeffs.ζₜ
    n₀², nₜ² = coeffs.n₀², coeffs.nₜ²
    h1₀, h2₀, h1ₜ, h2ₜ = coeffs.h1₀, coeffs.h2₀, coeffs.h1ₜ, coeffs.h2ₜ
    h1₀′, h2₀′, h1ₜ′, h2ₜ′ = coeffs.h1₀′, coeffs.h2₀′, coeffs.h1ₜ′, coeffs.h2ₜ′
    B1₀, B2₀, B3₀, B4₀ = coeffs.B1₀, coeffs.B2₀, coeffs.B3₀, coeffs.B4₀
    B1ₜ, B2ₜ, B3ₜ, B4ₜ = coeffs.B1ₜ, coeffs.B2ₜ, coeffs.B3ₜ, coeffs.B4ₜ
    A1, A2, A3, A4 = coeffs.A1, coeffs.A2, coeffs.A3, coeffs.A4

    n11, n22, d11, d22 = fresnelnd(C, W, ng²)
    dn11dC, dn22dC, dd11dC, dd22dC = dfresnelnddC(C, W, ng²)

    h1″ = -ζ₀*h1₀
    h2″ = -ζ₀*h2₀
    dzdC = 2*C*cbrtkoverαsq

    dB1dC = h1₀ + (C*h1₀′ + (K*h1″ + L*h1₀′)/n₀²)*dzdC
    dB2dC = h2₀ + (C*h2₀′ + (K*h2″ + L*h2₀′)/n₀²)*dzdC
    dB3dC = h1₀ + (C*h1₀′ + K*h1″)*dzdC
    dB4dC = h2₀ + (C*h2₀′ + K*h2″)*dzdC

    dA1dC = n11*dB1dC + dn11dC*B1₀ - 2(d11*h1₀′*dzdC + dd11dC*h1₀)
    dA2dC = n11*dB2dC + dn11dC*B2₀ - 2(d11*h2₀′*dzdC + dd11dC*h2₀)
    dA3dC = n22*dB3dC + dn22dC*B3₀ - 2(d22*h1₀′*dzdC + dd22dC*h1₀)
    dA4dC = n22*dB4dC + dn22dC*B4₀ - 2(d22*h2₀′*dzdC + dd22dC*h2₀)

    h1″ = -ζₜ*h1ₜ
    h2″ = -ζₜ*h2ₜ
    dzdC = 2*C*cbrtkoverαsq

    dB1dC = h1ₜ + (C*h1ₜ′ + (K*h1″ + L*h1ₜ′)/nₜ²)*dzdC
    dB2dC = h2ₜ + (C*h2ₜ′ + (K*h2″ + L*h2ₜ′)/nₜ²)*dzdC
    dB3dC = h1ₜ + (C*h1ₜ′ + K*h1″)*dzdC
    dB4dC = h2ₜ + (C*h2ₜ′ + K*h2″)*dzdC

    dd11dC = A2*dB1dC + dA2dC*B1ₜ - A1*dB2dC - dA1dC*B2ₜ
    dd22dC = A4*dB3dC + dA4dC*B3ₜ - A3*dB4dC - dA3dC*B4ₜ
    dn11dC = 2(A2*h1ₜ′*dzdC + dA2dC*h1ₜ - A1*h2ₜ′*dzdC - dA1dC*h2ₜ)
    dn22dC = 2(A4*h1ₜ′*dzdC + dA4dC*h1ₜ - A3*h2ₜ′*dzdC - dA3dC*h2ₜ)

    return dn11dC, dn22dC, dd11dC, dd22dC
end

"""
Return the Fresnel reflection matrix at the ground in terms of `n` and `d`.

See Morfitt Shellman 1976, Appendix C for the derivation.

See also: [`fresnelgroundcoeffs`](@ref)
"""
function fresnelnd(C, W, ng²)
    n11 = 1
    n22 = 1/W
    d11 = (C - W/ng²)/2
    d22 = (C/W - 1)/2

    return n11, n22, d11, d22
end

"""
Return the derivative of the Fresnel ground reflection terms `n` and `d` wrt cos(θ).

See also: [`fresnelgroundcoeffs`](@ref)
"""
function dfresnelnddC(C, W, ng²)
    W³ = W^3

    dn11dC = zero(ng²)
    dn22dC = -C/W³
    dd11dC = (1 - C/(W*ng²))/2
    dd22dC = (-C^2/W³ + 1/W)/2

    return dn11dC, dn22dC, dd11dC, dd22dC
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

See also: `mf_fctval.for`, `ndground`(@ref)
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
    groundcoeffs = _groundcoeffs(θ, ω, k, σ, ϵᵣ, reflectionheight, referenceheight)
    N11, N22, D11, D22 = ndground(groundcoeffs)

    k1 = N11 - X[1,1]*D11
    k2 = N22 - X[2,2]*D22
    f = k1*k2 - X[2,1]*X[1,2]*D11*D22
end

"""
Calculate modal excitation factors at reference height.

TODO: Support horizontal end and broadside exciters (antenna).

See Morfitt Shellman 1976 pg 29.
"""
function excitationfactors(θ)
    # Initialize
    S = sind(θ)

    dFdθ = ()

    B₁ = sqrt(S)^5/dFdθ
    B₂ = -B₁/S

    Ez = B₁*(1 + Rg11)^2*(1 - Rg22*R22)/(Rg11*D11)
    Ey = -B₁/S*R21*(1 + Rg11)*(1 + Rg22)/D12
    Ex = B₁/S*(1 + Rg11)^2*(1 - R22^2)/(Rg11*D11)  # TODO: Check this one

    return Ez, Ey, Ex
end

"""
Calculate modal height gain terms.
"""
function heightgain(coeffs::GroundCoefficients)
    C = cosd(θ)

    ζ = cbrt(k/α)^2*(C^2 + α*(height - reflectionheight))
    h1, h2, h1′, h2′ = modifiedhankel(ζ)

    fvertical(z) = (A2*h1ₜ - A1*h2ₜ)/-K
    fhorizontal(z) = (A4*h1ₜ - A3*h2ₜ)/(-W/K)
    hg(z) = -im*exp((z - d)/a)*(cbrt(2/(a*k))*(QC1*mh1p + GC1*mh2p) + (2/(a*k))*(QC1*mh1 + GC1*mh2))

    return fvertical, fhorizontal
end
