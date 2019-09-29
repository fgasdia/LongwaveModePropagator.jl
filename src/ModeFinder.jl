using SpecialFunctions
using LinearAlgebra
using StaticArrays
using PolynomialRoots
using DifferentialEquations

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

struct EigenAngle{T}
    θ::T
    cosθ::T
    sinθ::T
    cos²θ::T
    sin²θ::T

    function EigenAngle{T}(θ::T) where T <: Number
        C = cosd(θ)
        S = sind(θ)
        C² = C^2
        S² = 1 - C²
        new(θ, C, S, C², S²)
    end
end
EigenAngle(θ::T) where T <: Number = EigenAngle{T}(θ)

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
Computation of susceptibility `M` matrix as defined by Budden (1955)
[](10.1098/rspa.1955.0027).

# TODO: Add M up for each species
"""
function susceptibility(ω, z₀, z, spec::Constituent, bfield::BField)
    # Unpack
    B, l, m, n = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    l², m², n² = l^2, m^2, n^2

    # Constitutive relations (see Budden 1955, pg. 517)
    e, m, N, ν = spec.charge, spec.mass, spec.numberdensity, spec.collisionfrequency
    X = N(z)*e^2/(ϵ₀*m*ω^2)
    Y = e*B/(m*ω)
    Z = ν(z)/ω
    U = 1 - im*Z

    U² = U^2
    Y² = Y^2

    earthcurvature = 2(z₀ - z)/earthradius

    # In LWPC and Sheddy 1968 Fortran Program `earthcurvature` is not multiplied by capd
    # (the above line), even though it _is_ multiplied in MS 1976.
    # This seems to be supported by Pappert 1968
    M11 = U² - l²*Y²
    M21 = im*n*Y*U - l*m*Y²
    M31 = -im*m*Y*U - l*n*Y²
    M12 = -im*n*Y*U - l*m*Y²
    M22 = U² - m²*Y²
    M32 = im*l*Y*U - m*n*Y²
    M13 = im*m*Y*U - l*n*Y²
    M23 = -im*l*Y*U - m*n*Y²
    M33 =  U² - n²*Y²

    M = SMatrix{3,3}(M11, M21, M31,
                     M12, M22, M32,
                     M13, M23, M33)

    M *= -X/(U*(U² - Y²))

    M -= earthcurvature*I  # This correction only occurs once after adding species

    return M
end

"""
Compute matrix elements for solving differential of reflection matrix `R` wrt `z`.

See Budden 1955 second method.

```math
e′ = -iTe
```
"""
function smatrix(ea::EigenAngle, M)
    # Unpack
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = 1/C

    den = 1/(1 + M[3,3])

    # Temporary matrix elements T
    T11 = -S*M[3,1]*den
    T12 = S*M[3,2]*den
    # T13 = 0
    T14 = (C² + M[3,3])*den
    # T21 = 0
    # T22 = 0
    # T23 = 1
    # T24 = 0
    T31 = M[2,3]*M[3,1]*den - M[2,1]
    T32 = C² + M[2,2] - M[2,3]*M[3,2]*den
    # T33 = 0
    T34 = S*M[2,3]*den
    T41 = 1 + M[1,1] - M[1,3]*M[3,1]*den
    T42 = M[3,2]*M[1,3]*den - M[1,2]
    # T43 = 0
    T44 = -S*M[1,3]*den

    # Form the four 2x2 submatrices of `S`
    S11 = @SMatrix [T11+T44+T14*Cinv+C*T41 -T12*Cinv-T42;
                    -T31-T34*Cinv C+T32*Cinv]
    S12 = @SMatrix [-T11+T44+T14*Cinv-C*T41 -T12*Cinv-T42;
                    T31-T34*Cinv -C+T32*Cinv]
    S21 = @SMatrix [-T11+T44-T14*Cinv+C*T41 T12*Cinv-T42;
                    T31+T34*Cinv C-T32*Cinv]
    S22 = @SMatrix [T11+T44-T14*Cinv-C*T41 T12*Cinv-T42;
                    -T31+T34*Cinv -C-T32*Cinv]

    return S11, S12, S21, S22

    #== Equivalent to (but ~2x faster than):
    T = SMatrix{4,4}(T11, T21, T31, T41,
                     T12, T22, T32, T42,
                     T13, T23, T33, T43,
                     T14, T24, T34, T44)

    L = SMatrix{4,4}(C, 0, 0, 1,
                     0, -1, -C, 0,
                     -C, 0, 0, 1,
                     0, -1, C, 0)

    Linv = SMatrix{4,4}(Cinv, 0, -Cinv, 0,
                        0, -1, 0, -1,
                        0, -Cinv, 0, Cinv,
                        1, 0, 1, 0)

    S = Linv*T*L
    ==#
end

"""
Calculate the derivative of the reflection matrix `R` wrt height `z`.

Expects S to be a tuple of `(S11, S12, S21, S22)`.

The factor of `k` appears explicitly because Budden 1955 derives R′ wrt a height
variable ``s'' which includes k.
"""
function dRdz(R, S, k)
    return -im*k/2*(S[3] + S[4]*R - R*S[1] - R*S[2]*R)
end




#######################################################


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
function smatrix(ea::EigenAngle, M)
    # Unpack
    c, s, c² = ea.cosθ, ea.sinθ, ea.cos²θ

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

function dsmatrixdC(ea::EigenAngle, M)
    c, s, c² = ea.cosθ, ea.sinθ, ea.cos²θ

    ds = -c/s  # ds/dc
    dcs = s - c²/s

    oneplusM33 = 1 + M[3,3]

    T31 = M[2,3]*M[3,1]/oneplusM33 - M[2,1]
    T41 = 1 + M[1,1] - M[1,3]*M[3,1]/oneplusM33
    T42 = M[3,2]*M[1,3]/oneplusM33 - M[1,2]

    dB = @SMatrix [-2(T41+M[1,3]*ds/oneplusM33)   0;
                   0                              -2]
    dC = @SMatrix [-2(T41-M[3,1]*ds/oneplusM33)   0;
                   0                              -2]
    dD = @SMatrix [2*c*(T41-1+M[3,3]/oneplusM33)+dcs*(-M[3,1]+M[1,3])/oneplusM33    ds*M[3,2]/oneplusM33+T42;
                   -T31+ds*M[2,3]/oneplusM33                                        0]

    return dB, dC, dD
end

"""
Calculate angle from 315°.

Used to sort in [`sharplyboundedX!`](@ref).
"""
function anglefrom315(qval)
    angq = rad2deg(angle(qval))
    angq < 0 && (angq += 360)
    angq < 135 && (angq += 360)
    abs(angq - 315)
end

"""
Calculation of Booker quartic for solution of `X` for a sharply bounded ionosphere.

See also: [`sharplyboundedX!`](@ref)
"""
function bookerquartic(ea::EigenAngle, M)
    # Initialize
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

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

    bcoeffs = SVector{5}(b0, b1, b2, b3, b4)

    # TODO: (Once supported- MVector stays as MVector through ops) input bcoeffs to roots
    q = MVector{4}(roots([b0, b1, b2, b3, b4]))
    sort!(q, by=anglefrom315)

    return q, bcoeffs
end

"""
Common variables for [`sharplyboundedX`](@ref) and [`sharplyboundedXdXdC`](@ref).
"""
function _sharplybounded(ea::EigenAngle, M, q)
    # Initialize
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # Dispersion matrix elements for X
    D12 = M[1,2]
    D32 = M[3,2]
    D33 = C² + M[3,3]

    # For `q[1]`
    D11₁ = 1 + M[1,1] - q[1]^2
    D13₁ = M[1,3] + q[1]*S
    D31₁ = M[3,1] + q[1]*S
    denom₁ = D11₁*D33 - D13₁*D31₁
    P1 = (-D12*D33 + D13₁*D32)/denom₁
    T1 = q[1]*P1 - S*(-D11₁*D32 + D12*D31₁)/denom₁

    # For `q[2]`
    D11₂ = 1 + M[1,1] - q[2]^2
    D13₂ = M[1,3] + q[2]*S
    D31₂ = M[3,1] + q[2]*S
    denom₂ = D11₂*D33 - D13₂*D31₂
    P2 = (-D12*D33 + D13₂*D32)/denom₂
    T2 = q[2]*P2 - S*(-D11₂*D32 + D12*D31₂)/denom₂

    Δ = (T1*C + P1)*(C + q[2]) - (T2*C + P2)*(C + q[1])

    return (D12=D12, D32=D32, D33=D33,
            D11₁=D11₁, D13₁=D13₁, D31₁=D31₁, denom₁=denom₁, P1=P1, T1=T1,
            D11₂=D11₂, D13₂=D13₂, D31₂=D31₂, denom₂=denom₂, P2=P2, T2=T2,
            Δ=Δ)
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

See also: [`bookerquartic`](@ref)
"""
function sharplyboundedX!(X, D, ea::EigenAngle, q)
    # Unpack
    C = ea.cosθ
    P1, P2, T1, T2, Δ = D.P1, D.P2, D.T1, D.T2, D.Δ

    X[1,1] = T1*(C + q[2]) - T2*(C +q[1])
    X[2,1] = -(q[1] - q[2])
    X[1,2] = -(T1*P2 - T2*P1)
    X[2,2] = (T1*C + P1) - (T2*C + P2)

    X .*= 2/Δ

    return X
end

"""
There is likely an error in the LWPC code for this, `mf_initr.for` where `den` from `k=2` is
used for both `dq`.

TODO: Can probably do a matrix algebra version of this.
"""
function sharplyboundeddXdC!(dXdC, X, D, ea::EigenAngle, M, q, b)
    # Unpack
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ
    b0, b1, b2, b3, b4 = b

    D12, D32, D33 = D.D12, D.D32, D.D33
    D11₁, D13₁, D31₁ = D.D11₁, D.D13₁, D.D31₁
    D11₂, D13₂, D31₂ = D.D11₂, D.D13₂, D.D31₂
    denom₁, denom₂, P1, P2, T1, T2, Δ = D.denom₁, D.denom₂, D.P1, D.P2, D.T1, D.T2, D.Δ

    # Initialize
    dS = -C/S  # dS/dC

    # Booker quartic derivatives wrt C
    db3 = dS*(M[1,3] + M[3,1])
    db2 = -2C*(2 + M[1,1] + M[3,3])
    db1 = (dS/S)*b1 - 2S*C*(M[1,3] + M[3,1])
    db0 = 2C*((1 + M[1,1])*(C² + M[2,2] + C² + M[3,3]) - M[1,3]*M[3,1] - M[1,2]*M[2,1])

    # For `q[1]`
    dq1 = -(((db3*q[1] + db2)*q[1] + db1)*q[1] + db0) /
            (((4b4*q[1] + 3b3)*q[1] + 2b2)*q[1] + b1)

    dD11₁ = -2q[1]*dq1
    dD13₁ = q[1]*dS + dq1*S
    dD33 = 2C
    dD31₁ = dD13₁

    ddenom₁ = D11₁*dD33 + dD11₁*D33 - D13₁*dD31₁ - dD13₁*D31₁
    dP1 = -(ddenom₁*P1 + (D12*dD33 - dD13₁*D32))/denom₁
    dT1 = q[1]*dP1 + dq1*P1 +
        (D11₁*D32 - D12*D31₁)*(dS - S*ddenom₁/denom₁)/denom₁ +
        S*(dD11₁*D32 - D12*dD31₁)/denom₁  # XXX: LWPC `mf_initr.for` likely has a bug, not using the correct denom

    # For `q[2]`
    dq2 = -(((db3*q[2] + db2)*q[2] + db1)*q[2] + db0) /
            (((4b4*q[2] + 3b3)*q[2] + 2b2)*q[2] + b1)

    dD11₂ = -2q[2]*dq2
    dD13₂ = q[2]*dS + dq2*S
    dD31₂ = dD13₂

    ddenom₂ = D11₂*dD33 + dD11₂*D33 - D13₂*dD31₂ - dD13₂*D31₂
    dP2 = -(ddenom₂*P2 + (D12*dD33 - dD13₂*D32))/denom₂
    dT2 = q[2]*dP2 + dq2*P2 +
        (D11₂*D32 - D12*D31₂)*(dS - S*ddenom₂/denom₂)/denom₂ +
        S*(dD11₂*D32 - D12*dD31₂)/denom₂

    dΔ = (T1*C + P1)*(1 + dq2) +
         (T1 + dT1*C + dP1)*(C + q[2]) -
         (T2*C + P2)*(1 + dq1) -
         (T2 + dT2*C + dP2)*(C + q[1])

    factor1 = 2/Δ
    factor2 = dΔ/Δ

    dXdC[1,1] = -X[1,1]*factor2 +
        factor1*(T1*(1 + dq2) + dT1*(C + q[2]) - T2*(1 + dq1) - dT2*(C + q[1]))
    dXdC[2,1] = -X[2,1]*factor2 -
        factor1*(dq1 - dq2)
    dXdC[1,2] = -X[1,2]*factor2 -
        factor1*(T1*dP2 + dT1*P2 - T2*dP1 - dT2*P1)
    dXdC[2,2] = -X[2,2]*factor2 +
        factor1*((T1 + dT1*C + dP1) - (T2 + dT2*C + dP2))

    return dXdC
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
function dXdCdh(X, dX, B, dB, C, dC, D, dD, k)
    -im*k/2*(B*dX + dB*X + X*dC + dX*C + X*D*dX + X*dD*X + dX*D*X)
end

"""
"""
function Xderivative(X, params, height)
    M = mmatrix!(params.M, params.ω, params.referenceheight, height,
                 params.species, params.bfield)
    A, B, C, D = smatrix(params.ea, M)
    dXdh(X, A, B, C, D, params.k)
end

function XdXdCderivative(XdX, params, height)
    M = mmatrix!(params.M, params.ω, params.referenceheight, height,
                 params.species, params.bfield)
    A, B, C, D = smatrix(params.ea, M)
    dB, dC, dD = dsmatrixdC(params.ea, M)
    X = @view XdX[1:2,:]
    dX = @view XdX[3:4,:]
    vcat(dXdh(X, A, B, C, D, params.k), dXdCdh(X, dX, B, dB, C, dC, D, dD, params.k))
end

"""
Full wave solution for `X` through ionosphere.

Integrates ``dX/dh`` from `topheight` to `bottomheight`.
"""
function integratethroughionosphere(ea::EigenAngle, source::AbstractSource, fromheight, toheight, referenceheight,
                                    species, bfield::BField)

    # Unpack
    ω, k = source.ω, source.k  # k in km

    # Initialize
    M = @MMatrix zeros(ComplexF64, 3, 3)
    X = @MMatrix zeros(ComplexF64, 2, 2)

    # Initial condition for integration
    # `fromheight` should be topheight
    mmatrix!(M, ω, referenceheight, fromheight, species, bfield)
    q, bcoeffs = bookerquartic(ea, M)
    D = _sharplybounded(ea, M, q)
    initialX = sharplyboundedX!(X, D, ea, q)

    params = (M=M, ea=ea, referenceheight=referenceheight, ω=ω, k=k, species=species,
              bfield=bfield)

    heightbounds = (fromheight, toheight)
    prob = ODEProblem(Xderivative, initialX, heightbounds, params)
    # sol = solve(prob, BS3(), reltol=1e-3, dtmax=1., save_everystep=false, save_start=false)  # Faster, less accurate
    sol = solve(prob, Tsit5(), reltol=1e-4, dtmax=1., save_everystep=false, save_start=false)  # Decent accuracy for the speed

    return sol
end

function integratethroughionosphere_dC(ea::EigenAngle, source::AbstractSource, fromheight, toheight, referenceheight,
                                       species, bfield::BField)
    # Unpack
    ω, k = source.ω, source.k  # k in km

    # Initialize
    M = @MMatrix zeros(ComplexF64, 3, 3)
    X = @MMatrix zeros(ComplexF64, 2, 2)
    dX = similar(X)

    mmatrix!(M, ω, referenceheight, fromheight, species, bfield)
    q, bcoeffs = bookerquartic(ea, M)
    D = _sharplybounded(ea, M, q)
    initialX = sharplyboundedX!(X, D, ea, q)
    initialdX = sharplyboundeddXdC!(dX, X, D, ea, M, q, bcoeffs)

    params = (M=M, ea=ea, referenceheight=referenceheight, ω=ω, k=k, species=species,
              bfield=bfield)

    heightbounds = (fromheight, toheight)
    prob = ODEProblem(XdXdCderivative, vcat(initialX, initialdX), heightbounds, params)

    sol = solve(prob, AutoTsit5(Rosenbrock23()), save_everystep=false, save_start=false)

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

Boundary conditions are set at the bottom so we can solve for ``(-A/C)``, ``(-B/C)``,
``(-Q/C)``, and ``(-G/C)``. Then we use those values to calculate the ``E`` fields at
`toheight`.
"""
function integratethroughfreespace!(X, ea::EigenAngle, k, fromheight, toheight, referenceheight)
    # Unpack
    C, C² = ea.cosθ, ea.cos²θ

    # Initialize
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2k)
    cbrtkoverαsq = cbrt(k/α)^2  # (k/α)^(2/3)

    # Complex "heights"
    ζ₀ = cbrtkoverαsq*(C² + α*(fromheight - referenceheight))
    n₀² = 1 + α*(fromheight - referenceheight)

    # Computation of height-gain coefficients for two conditions on the upgoing wave from
    # `fromheight`; E∥ = 1, Ey = 0 and E∥ = 0, Ey = 1
    h1₀, h2₀, h1′₀, h2′₀ = modifiedhankel(ζ₀)

    # Horizontal polarization, BC 1 and 2
    # NOTE: There is no real advantage to having `a` be even an SMatrix
    aₕ11 = h1₀
    aₕ12 = h2₀
    aₕ21 = C*h1₀ + K*h1′₀
    aₕ22 = C*h2₀ + K*h2′₀
    Δₕ = aₕ11*aₕ22 - aₕ21*aₕ12

    AC₁ = X[2,1]*aₕ22/Δₕ
    BC₁ = -X[2,1]*aₕ21/Δₕ
    AC₂ = (X[2,2]*aₕ22 - 2aₕ12)/Δₕ
    BC₂ = (2aₕ11 - X[2,2]*aₕ21)/Δₕ

    # Vertical polarization, boundary conditions 1 and 2
    aᵥ11 = h1₀
    aᵥ12 = h2₀
    aᵥ21 = C*h1₀ + (K*h1′₀ + L*h1₀)/n₀²
    aᵥ22 = C*h2₀ + (K*h2′₀ + L*h2₀)/n₀²
    Δᵥ = aᵥ11*aᵥ22 - aᵥ21*aᵥ12

    QC₁ = (X[1,1]*aᵥ22 - 2aᵥ12)/Δᵥ
    GC₁ = (2aᵥ11 - X[1,1]*aᵥ21)/Δᵥ
    QC₂ = X[1,2]*aᵥ22/Δᵥ
    GC₂ = -X[1,2]*aᵥ21/Δᵥ

    # Computation of upgoing fields E∥ and Ey at `toheight` for the two conditions above.
    ζₜ = cbrtkoverαsq*(C² + α*(toheight - referenceheight))
    nₜ² = 1 + α*(toheight - referenceheight)

    h1ₜ, h2ₜ, h1′ₜ, h2′ₜ = modifiedhankel(ζₜ)

    aₕₜ21 = C*h1ₜ + K*h1′ₜ
    aₕₜ22 = C*h2ₜ + K*h2′ₜ

    # Calculate parallel (p) and y fields
    Eyt₁ = (AC₁*aₕₜ21 + BC₁*aₕₜ22)/2
    Eyt₂ = (AC₂*aₕₜ21 + BC₂*aₕₜ22)/2

    aᵥₜ21 = C*h1ₜ + (K*h1′ₜ + L*h1ₜ)/nₜ²
    aᵥₜ22 = C*h2ₜ + (K*h2′ₜ + L*h2ₜ)/nₜ²

    Ept₁ = (QC₁*aᵥₜ21 + GC₁*aᵥₜ22)/2
    Ept₂ = (QC₂*aᵥₜ21 + GC₂*aᵥₜ22)/2

    # Reflection matrix at the `toheight` level
    W = Ept₁*Eyt₂ - Ept₂*Eyt₁
    V₁h = AC₁*h1ₜ + BC₁*h2ₜ
    V₂h = AC₂*h1ₜ + BC₂*h2ₜ
    V₁v = QC₁*h1ₜ + GC₁*h2ₜ
    V₂v = QC₂*h1ₜ + GC₂*h2ₜ

    # Mutate to `X` at the `toheight` level
    X[1,1] = V₁v*Eyt₂ - V₂v*Eyt₁
    X[2,1] = V₁h*Eyt₂ - V₂h*Eyt₁
    X[1,2] = V₂v*Ept₁ - V₁v*Ept₂
    X[2,2] = V₂h*Ept₁ - V₁h*Ept₂
    X ./= W

    # X = @SMatrix [(V₁v*Eyt₂ - V₂v*Eyt₁)/W   (V₂v*Ept₁ - V₁v*Ept₂)/W;
    #               (V₁h*Eyt₂ - V₂h*Eyt₁)/W   (V₂h*Ept₁ - V₁h*Ept₂)/W]

    return X
end

"""
NOTE: This function cannot be separated from the X update because we need both `X from` and
`X to` to solve for `dXdC from` and `dXdC to` but we update `X` in place.
"""
function integratethroughfreespace_XdXdC!(X, dXdC, ea::EigenAngle, k, fromheight, toheight, referenceheight)
    # Unpack
    C, C² = ea.cosθ, ea.cos²θ

    # Initialize
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2k)
    cbrtkoverαsq = cbrt(k/α)^2  # (k/α)^(2/3)

    # Complex "heights"
    ζ₀ = cbrtkoverαsq*(C² + α*(fromheight - referenceheight))
    dζdC = 2C*cbrtkoverαsq  # irrespective of height
    n₀² = 1 + α*(fromheight - referenceheight)

    # Computation of height-gain coefficients for two conditions on the upgoing wave from
    # `fromheight`; E∥ = 1, Ey = 0 and E∥ = 0, Ey = 1
    h1₀, h2₀, h1′₀, h2′₀ = modifiedhankel(ζ₀)
    h1″₀ = -ζ₀*h1₀
    h2″₀ = -ζ₀*h2₀

    # Horizontal polarization, BC 1 and 2
    # NOTE: There is no real advantage to having `a` be even an SMatrix
    aₕ11 = h1₀
    aₕ12 = h2₀
    aₕ21 = C*h1₀ + K*h1′₀
    aₕ22 = C*h2₀ + K*h2′₀
    Δₕ = aₕ11*aₕ22 - aₕ21*aₕ12

    daₕ11dC = h1′₀*dζdC
    daₕ12dC = h2′₀*dζdC
    daₕ21dC = h1₀ + (C*h1′₀ + K*h1″₀)*dζdC
    daₕ22dC = h2₀ + (C*h2′₀ + K*h2″₀)*dζdC
    dΔₕdC = aₕ11*daₕ22dC + daₕ11dC*aₕ22 - aₕ21*daₕ12dC - daₕ21dC*aₕ12

    AC₁ = X[2,1]*aₕ22/Δₕ
    BC₁ = -X[2,1]*aₕ21/Δₕ
    AC₂ = (X[2,2]*aₕ22 - 2aₕ12)/Δₕ
    BC₂ = (2aₕ11 - X[2,2]*aₕ21)/Δₕ

    tmp1 = dΔₕdC/Δₕ
    dAC₁dC = (X[2,1]*daₕ22dC + dXdC[2,1]*aₕ22)/Δₕ - AC₁*tmp1
    dBC₁dC = (-X[2,1]*daₕ21dC - dXdC[2,1]*aₕ21)/Δₕ - BC₁*tmp1
    dAC₂dC = (X[2,2]*daₕ22dC + dXdC[2,2]*aₕ22 - 2daₕ12dC)/Δₕ - AC₂*tmp1
    dBC₂dC = (2daₕ11dC - X[2,2]*daₕ21dC - dXdC[2,2]*aₕ21)/Δₕ - BC₂*tmp1

    # Vertical polarization, boundary conditions 1 and 2
    aᵥ11 = h1₀
    aᵥ12 = h2₀
    aᵥ21 = C*h1₀ + (K*h1′₀ + L*h1₀)/n₀²
    aᵥ22 = C*h2₀ + (K*h2′₀ + L*h2₀)/n₀²
    Δᵥ = aᵥ11*aᵥ22 - aᵥ21*aᵥ12

    daᵥ11dC = h1′₀*dζdC
    daᵥ12dC = h2′₀*dζdC
    daᵥ21dC = h1₀ + (C*h1′₀ + (K*h1″₀ + L*h1′₀)/n₀²)*dζdC
    daᵥ22dC = h2₀ + (C*h2′₀ + (K*h2″₀ + L*h2′₀)/n₀²)*dζdC
    dΔᵥdC = aᵥ11*daᵥ22dC + daᵥ11dC*aᵥ22 - aᵥ21*daᵥ12dC - daᵥ21dC*aᵥ12

    QC₁ = (X[1,1]*aᵥ22 - 2aᵥ12)/Δᵥ
    GC₁ = (2aᵥ11 - X[1,1]*aᵥ21)/Δᵥ
    QC₂ = X[1,2]*aᵥ22/Δᵥ
    GC₂ = -X[1,2]*aᵥ21/Δᵥ

    tmp2 = dΔᵥdC/Δᵥ
    dQC₁dC = (X[1,1]*daᵥ22dC + dXdC[1,1]*aᵥ22 - 2daᵥ12dC)/Δᵥ - QC₁*tmp2
    dGC₁dC = (2daᵥ11dC - X[1,1]*daᵥ21dC - dXdC[1,1]*aᵥ21)/Δᵥ - GC₁*tmp2
    dQC₂dC = (X[1,2]*daᵥ22dC + dXdC[1,2]*aᵥ22)/Δᵥ - QC₂*tmp2
    dGC₂dC = (-X[1,2]*daᵥ21dC - dXdC[1,2]*aᵥ21)/Δᵥ - GC₂*tmp2

    # Computation of upgoing fields E∥ and Ey at `toheight` for the two conditions above.
    ζₜ = cbrtkoverαsq*(C² + α*(toheight - referenceheight))
    nₜ² = 1 + α*(toheight - referenceheight)

    h1ₜ, h2ₜ, h1′ₜ, h2′ₜ = modifiedhankel(ζₜ)
    h1″ₜ = -ζₜ*h1ₜ
    h2″ₜ = -ζₜ*h2ₜ

    aₕₜ21 = C*h1ₜ + K*h1′ₜ
    aₕₜ22 = C*h2ₜ + K*h2′ₜ

    daₕₜ21dC = h1ₜ + (C*h1′ₜ + K*h1″ₜ)*dζdC
    daₕₜ22dC = h2ₜ + (C*h2′ₜ + K*h2″ₜ)*dζdC

    # Calculate parallel (p) and y fields
    Eyt₁ = (AC₁*aₕₜ21 + BC₁*aₕₜ22)/2
    Eyt₂ = (AC₂*aₕₜ21 + BC₂*aₕₜ22)/2

    dEyt₁dC = (AC₁*daₕₜ21dC + dAC₁dC*aₕₜ21 + BC₁*daₕₜ22dC + dBC₁dC*aₕₜ22)/2
    dEyt₂dC = (AC₂*daₕₜ21dC + dAC₂dC*aₕₜ21 + BC₂*daₕₜ22dC + dBC₂dC*aₕₜ22)/2

    aᵥₜ21 = C*h1ₜ + (K*h1′ₜ + L*h1ₜ)/nₜ²
    aᵥₜ22 = C*h2ₜ + (K*h2′ₜ + L*h2ₜ)/nₜ²

    daᵥₜ21dC = h1ₜ + (C*h1′ₜ + (K*h1″ₜ + L*h1′ₜ)/nₜ²)*dζdC
    daᵥₜ22dC = h2ₜ + (C*h2′ₜ + (K*h2″ₜ + L*h2′ₜ)/nₜ²)*dζdC

    Ept₁ = (QC₁*aᵥₜ21 + GC₁*aᵥₜ22)/2
    Ept₂ = (QC₂*aᵥₜ21 + GC₂*aᵥₜ22)/2

    dEpt₁dC = (QC₁*daᵥₜ21dC + dQC₁dC*aᵥₜ21 + GC₁*daᵥₜ22dC + dGC₁dC*aᵥₜ22)/2
    dEpt₂dC = (QC₂*daᵥₜ21dC + dQC₂dC*aᵥₜ21 + GC₂*daᵥₜ22dC + dGC₂dC*aᵥₜ22)/2

    # Reflection matrix at the `toheight` level
    W = Ept₁*Eyt₂ - Ept₂*Eyt₁
    V₁h = AC₁*h1ₜ + BC₁*h2ₜ
    V₂h = AC₂*h1ₜ + BC₂*h2ₜ
    V₁v = QC₁*h1ₜ + GC₁*h2ₜ
    V₂v = QC₂*h1ₜ + GC₂*h2ₜ

    dWdC = Ept₁*dEyt₂dC + dEpt₁dC*Eyt₂ - Ept₂*dEyt₁dC - dEpt₂dC*Eyt₁
    dV₁hdC = (AC₁*h1′ₜ + BC₁*h2′ₜ)*dζdC + dAC₁dC*h1ₜ + dBC₁dC*h2ₜ
    dV₂hdC = (AC₂*h1′ₜ + BC₂*h2′ₜ)*dζdC + dAC₂dC*h1ₜ + dBC₂dC*h2ₜ
    dV₁vdC = (QC₁*h1′ₜ + GC₁*h2′ₜ)*dζdC + dQC₁dC*h1ₜ + dGC₁dC*h2ₜ
    dV₂vdC = (QC₂*h1′ₜ + GC₂*h2′ₜ)*dζdC + dQC₂dC*h1ₜ + dGC₂dC*h2ₜ

    # Mutate to `X` and `dXdC` at the `toheight`
    X[1,1] = V₁v*Eyt₂ - V₂v*Eyt₁
    X[2,1] = V₁h*Eyt₂ - V₂h*Eyt₁
    X[1,2] = V₂v*Ept₁ - V₁v*Ept₂
    X[2,2] = V₂h*Ept₁ - V₁h*Ept₂
    X ./= W

    tmp3 = dWdC/W
    dXdC[1,1] = (V₁v*dEyt₂dC + dV₁vdC*Eyt₂ - V₂v*dEyt₁dC - dV₂vdC*Eyt₁)/W - X[1,1]*tmp3
    dXdC[2,1] = (V₁h*dEyt₂dC + dV₁hdC*Eyt₂ - V₂h*dEyt₁dC - dV₂hdC*Eyt₁)/W - X[2,1]*tmp3
    dXdC[1,2] = (V₂v*dEpt₁dC + dV₂vdC*Ept₁ - V₁v*dEpt₂dC - dV₁vdC*Ept₂)/W - X[1,2]*tmp3
    dXdC[2,2] = (V₂h*dEpt₁dC + dV₂hdC*Ept₁ - V₁h*dEpt₂dC - dV₁hdC*Ept₂)/W - X[2,2]*tmp3

    return X, dXdC
end

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
function _groundcoeffs(ea::EigenAngle, source::AbstractSource, σ, ϵᵣ, toheight, referenceheight)
    # Unpack
    ω, k = source.ω, source.k
    C, C², S² = ea.cosθ, ea.cos²θ, ea.sin²θ

    # Initialize
    α = 2/earthradius
    K = im*cbrt(α/k)
    L = im*(α/2k)
    cbrtkoverαsq = cbrt(k/α)^2

    # At the "from" level, ``z = 0`` (ground)
    ng² = complex(ϵᵣ, -σ/(ω*ϵ₀))
    W = sqrt(ng² - S²)

    n11, d11, n22, d22 = fresnelnd(C, W, ng²)

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

    A1 = n11*B1₀ - 2d11*h1₀
    A2 = n11*B2₀ - 2d11*h2₀
    A3 = n22*B3₀ - 2d22*h1₀
    A4 = n22*B4₀ - 2d22*h2₀

    B1ₜ = C*h1ₜ + (K*h1ₜ′ + L*h1ₜ)/nₜ²
    B2ₜ = C*h2ₜ + (K*h2ₜ′ + L*h2ₜ)/nₜ²
    B3ₜ = C*h1ₜ + K*h1ₜ′
    B4ₜ = C*h2ₜ + K*h2ₜ′

    return (W=W, ng²=ng², K=K, L=L, cbrtkoverαsq=cbrtkoverαsq, ζ₀=ζ₀, ζₜ=ζₜ, n₀²=n₀², nₜ²=nₜ²,
            B1₀=B1₀, B2₀=B2₀, B3₀=B3₀, B4₀=B4₀, B1ₜ=B1ₜ, B2ₜ=B2ₜ, B3ₜ=B3ₜ, B4ₜ=B4ₜ,
            A1=A1, A2=A2, A3=A3, A4=A4,
            h1₀=h1₀, h2₀=h2₀, h1₀′=h1₀′, h2₀′=h2₀′, h1ₜ=h1ₜ, h2ₜ=h2ₜ, h1ₜ′=h1ₜ′, h2ₜ′=h2ₜ′)
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
function ndground(coeffs)
    A1, A2, A3, A4 = coeffs.A1, coeffs.A2, coeffs.A3, coeffs.A4
    B1ₜ, B2ₜ, B3ₜ, B4ₜ = coeffs.B1ₜ, coeffs.B2ₜ, coeffs.B3ₜ, coeffs.B4ₜ
    h1ₜ, h2ₜ = coeffs.h1ₜ, coeffs.h2ₜ

    N11 = 2(A2*h1ₜ - A1*h2ₜ)
    N22 = 2(A4*h1ₜ - A3*h2ₜ)
    D11 = A2*B1ₜ - A1*B2ₜ
    D22 = A4*B3ₜ - A3*B4ₜ

    return N11, D11, N22, D22
end

"""
Return the derivative of the ground reflection coefficient terms `n` and `d` wrt θ at the
reflection height.

See also: [`ndground`](@ref)
"""
function dndgrounddC(ea::EigenAngle, coeffs)
    # Unpack
    C, C² = ea.cosθ, ea.cos²θ
    W, ng² = coeffs.W, coeffs.ng²
    K, L, cbrtkoverαsq = coeffs.K, coeffs.L, coeffs.cbrtkoverαsq
    ζ₀, ζₜ = coeffs.ζ₀, coeffs.ζₜ
    n₀², nₜ² = coeffs.n₀², coeffs.nₜ²
    h1₀, h2₀, h1ₜ, h2ₜ = coeffs.h1₀, coeffs.h2₀, coeffs.h1ₜ, coeffs.h2ₜ
    h1₀′, h2₀′, h1ₜ′, h2ₜ′ = coeffs.h1₀′, coeffs.h2₀′, coeffs.h1ₜ′, coeffs.h2ₜ′
    B1₀, B2₀, B3₀, B4₀ = coeffs.B1₀, coeffs.B2₀, coeffs.B3₀, coeffs.B4₀
    B1ₜ, B2ₜ, B3ₜ, B4ₜ = coeffs.B1ₜ, coeffs.B2ₜ, coeffs.B3ₜ, coeffs.B4ₜ
    A1, A2, A3, A4 = coeffs.A1, coeffs.A2, coeffs.A3, coeffs.A4

    n11, d11, n22, d22 = fresnelnd(C, W, ng²)
    dn11dC, dd11dC, dn22dC, dd22dC = dfresnelnddC(C, C², W, ng²)

    h1″ = -ζ₀*h1₀
    h2″ = -ζ₀*h2₀
    dζdC = 2C*cbrtkoverαsq

    dB1dC = h1₀ + (C*h1₀′ + (K*h1″ + L*h1₀′)/n₀²)*dζdC
    dB2dC = h2₀ + (C*h2₀′ + (K*h2″ + L*h2₀′)/n₀²)*dζdC
    dB3dC = h1₀ + (C*h1₀′ + K*h1″)*dζdC
    dB4dC = h2₀ + (C*h2₀′ + K*h2″)*dζdC

    dA1dC = n11*dB1dC + dn11dC*B1₀ - 2(d11*h1₀′*dζdC + dd11dC*h1₀)
    dA2dC = n11*dB2dC + dn11dC*B2₀ - 2(d11*h2₀′*dζdC + dd11dC*h2₀)
    dA3dC = n22*dB3dC + dn22dC*B3₀ - 2(d22*h1₀′*dζdC + dd22dC*h1₀)
    dA4dC = n22*dB4dC + dn22dC*B4₀ - 2(d22*h2₀′*dζdC + dd22dC*h2₀)

    h1″ = -ζₜ*h1ₜ
    h2″ = -ζₜ*h2ₜ

    dB1dC = h1ₜ + (C*h1ₜ′ + (K*h1″ + L*h1ₜ′)/nₜ²)*dζdC
    dB2dC = h2ₜ + (C*h2ₜ′ + (K*h2″ + L*h2ₜ′)/nₜ²)*dζdC
    dB3dC = h1ₜ + (C*h1ₜ′ + K*h1″)*dζdC
    dB4dC = h2ₜ + (C*h2ₜ′ + K*h2″)*dζdC

    dd11dC = A2*dB1dC + dA2dC*B1ₜ - A1*dB2dC - dA1dC*B2ₜ
    dd22dC = A4*dB3dC + dA4dC*B3ₜ - A3*dB4dC - dA3dC*B4ₜ
    dn11dC = 2(A2*h1ₜ′*dζdC + dA2dC*h1ₜ - A1*h2ₜ′*dζdC - dA1dC*h2ₜ)
    dn22dC = 2(A4*h1ₜ′*dζdC + dA4dC*h1ₜ - A3*h2ₜ′*dζdC - dA3dC*h2ₜ)

    return dn11dC, dd11dC, dn22dC, dd22dC
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

    return n11, d11, n22, d22
end

"""
Return the derivative of the Fresnel ground reflection terms `n` and `d` wrt cos(θ).

See also: [`fresnelgroundcoeffs`](@ref)
"""
function dfresnelnddC(C, C², W, ng²)
    W³ = W^3

    dn11dC = zero(ng²)
    dn22dC = -C/W³
    dd11dC = (1 - C/(W*ng²))/2
    dd22dC = (-C²/W³ + 1/W)/2

    return dn11dC, dd11dC, dn22dC, dd22dC
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
function modifiedmodalfunction(ea::EigenAngle, source::AbstractSource, σ, ϵᵣ,
                               bottomheight, topheight, reflectionheight, referenceheight,
                               species,
                               bfield::BField)
    # Unpack
    k = source.k

    # Calculate `X` at `bottomheight`
    Xsol = integratethroughionosphere(ea, source, topheight, bottomheight, referenceheight,
                                      species, bfield)
    X = Xsol[end]

    # Refer `X` from `bottomheight` to `referenceheight`
    integratethroughfreespace!(X, ea, k, bottomheight, reflectionheight, referenceheight)

    # Calculate ground reflection matrix as `N` and `D` referred to `referenceheight`
    groundcoeffs = _groundcoeffs(ea, source, σ, ϵᵣ, reflectionheight, referenceheight)
    N11, D11, N22, D22 = ndground(groundcoeffs)

    k1 = N11 - X[1,1]*D11
    k2 = N22 - X[2,2]*D22
    f = k1*k2 - X[2,1]*X[1,2]*D11*D22
end

function modifiedmodalfunction_XdXdθ(ea::EigenAngle, source::AbstractSource, σ, ϵᵣ,
                                 bottomheight, topheight, reflectionheight, referenceheight,
                                 species,
                                 bfield::BField)

    # Unpack
    k = source.k
    S = ea.sinθ

    # Calculate `X` and `dXdC` at `bottomheight`
    sol = integratethroughionosphere_dC(ea, source,
                                        topheight, bottomheight, referenceheight,
                                        species, bfield)
    X = @view sol[end][1:2,:]
    dX = @view sol[end][3:4,:]

    # Refer `X` and `dXdC` from `bottomheight` to `referenceheight`
    integratethroughfreespace_XdXdC!(X, dX, ea, k,
                                    bottomheight, reflectionheight, referenceheight)

    # Calculate ground reflection matrix as `N` and `D` referred to     referenceheight
    groundcoeffs = _groundcoeffs(ea, source, σ, ϵᵣ, reflectionheight, referenceheight)
    N11, D11, N22, D22 = ndground(groundcoeffs)
    dND = SMatrix{2,2}(dndgrounddC(ea, groundcoeffs))

    D11D22 = D11*D22
    X21X12 = X[2,1]*X[1,2]
    k1 = N11 - X[1,1]*D11  # = dfdN11
    k2 = N22 - X[2,2]*D22  # = dfdN22
    f = k1*k2 - X21X12*D11D22

    dfdND = @SMatrix [k2                        k1;
                      -X[1,1]*k2-X21X12*D22     -X[2,2]*k1-X21X12*D11]

    dfdX = @SMatrix [-D11*k2           -X[2,1]*D11D22;
                     -X[1,2]*D11D22    -D22*k1]

    df = sum(dfdX .* dX + dfdND .* dND)

    dCdθ = deg2rad(-S)
    dfdθ = df*dCdθ

    return f, dfdθ
end

R2X(R, C) = (R + I)/C
X2R(X, C) = C*X - I

function modeparameters(ea::EigenAngle, X)
    # Unpack
    C, C², S = ea.cosθ, ea.cos²θ, ea.sinθ

    R = X2R(X, C)

    D11 = C*N11ₘ - D11ₘ
    D22 = C*N22ₘ - D22ₘ

    Rg = Diagonal(SVector{2}(D11ₘ/D11, D22ₘ/D22))
    Rgp1 = C*Diagonal(SVector{2}(N11ₘ/D11, N22ₘ/D22))

    dfdθ = rad2deg(dfdθₘ*C²/(D11*D22))
    factor = sqrt(S)/dfdθ

    # Eigen angle referred to ground level
    K = sqrt(1 - 2h/earthradius)
    stp = S/K
    tp = rad2deg(-log(sqrt(1 - stp^2) + im*stp))

    N = I - Rg*R

    # Excitation terms

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
function heightgain(ea::EigenAngle, groundcoeffs)
    C = ea.cosθ

    ζ = cbrt(k/α)^2*(C^2 + α*(height - reflectionheight))
    h1, h2, h1′, h2′ = modifiedhankel(ζ)

    fvertical(z) = (A2*h1ₜ - A1*h2ₜ)/-K
    fhorizontal(z) = (A4*h1ₜ - A3*h2ₜ)/(-W/K)
    hg(z) = -im*exp((z - d)/a)*(cbrt(2/(a*k))*(QC1*mh1p + GC1*mh2p) + (2/(a*k))*(QC1*mh1 + GC1*mh2))

    return fvertical, fhorizontal
end
