using SpecialFunctions
using LinearAlgebra
using StaticArrays
using OrdinaryDiffEq

using PolynomialRoots: roots!

using ModifiedHankelFunctionsOfOrderOneThird

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

`z` is "current" height
`z₀` is reference height for earth curvature where the index of refraction ≈ 1

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

    M -= earthcurvature*I  # This correction only occurs once after adding species?

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
    m31d = M[3,1]*den
    m32d = M[3,2]*den

    T11 = -S*m31d
    T12 = S*m32d
    # T13 = 0
    T14 = (C² + M[3,3])*den
    # T21 = 0
    # T22 = 0
    # T23 = 1
    # T24 = 0
    T31 = M[2,3]*m31d - M[2,1]
    T32 = C² + M[2,2] - M[2,3]*m32d
    # T33 = 0
    T34 = S*M[2,3]*den
    T41 = 1 + M[1,1] - M[1,3]*m31d
    T42 = M[1,3]*m32d - M[1,2]
    # T43 = 0
    T44 = -S*M[1,3]*den

    # S matrix, based on Sheddy et al., 1968, A Fortran Program...
    t12Cinv = T12*Cinv
    t14Cinv = T14*Cinv
    t32Cinv = T32*Cinv
    t34Cinv = T34*Cinv
    t41C = C*T41

    s11a = T11 + T44
    d11a = T11 - T44
    s11b = t14Cinv + t41C
    d11b = t14Cinv - t41C
    s12 = t12Cinv + T42
    d12 = t12Cinv - T42
    s21 = T31 + t34Cinv
    d21 = T31 - t34Cinv
    s22 = C + t32Cinv
    d22 = C - t32Cinv

    # Form the four 2x2 submatrices of `S`
    S11 = @SMatrix [s11a+s11b -s12;
                    -s21 s22]
    S12 = @SMatrix [-d11a+d11b -s12;
                    d21 -d22]
    S21 = @SMatrix [-d11a-d11b d12;
                    s21 d22]
    S22 = @SMatrix [s11a-s11b d12;
                    -d21 -s22]

    return S11, S12, S21, S22

    #== Equivalent to (but faster than):
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
Calculation of Booker quartic for solution of `R` for a sharply bounded ionosphere.

Based on Sheddy 1968 A General Analytic Solution for Reflection from a Sharply...

See also: [`sharplyboundedX!`](@ref)
"""
function bookerquartic(ea::EigenAngle, M)
    # Initialize
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # Booker quartic coefficients
    B4 = 1 + M[3,3]
    B3 = S*(M[1,3] + M[3,1])
    B2 = -(C² + M[3,3])*(1 + M[1,1]) + M[1,3]*M[3,1] -
         (1 + M[3,3])*(C² + M[2,2]) + M[2,3]*M[3,2]
    B1 = S*(M[1,2]*M[2,3] + M[2,1]*M[3,2] -
         (C² + M[2,2])*(M[1,3] + M[3,1]))
    B0 = (1 + M[1,1])*(C² + M[2,2])*(C² + M[3,3]) +
         M[1,2]*M[2,3]*M[3,1] + M[1,3]*M[2,1]*M[3,2] -
         M[1,3]*(C² + M[2,2])*M[3,1] -
         (1 + M[1,1])*M[2,3]*M[3,2] -
         M[1,2]*M[2,1]*(C² + M[3,3])

    q = @MVector zeros(complex(typeof(S)), 4)  # algorithm requires complex type
    B = MVector{5}(B0, B1, B2, B3, B4)

    roots!(q, B, NaN, 4, false)

    return q
end

"""
Calculate angle from 315°.

Used to sort in [`sharplyboundedX!`](@ref).
"""
function anglefrom315(qval)
    angq = rad2deg(angle(qval))  # rad2deg adds a few ns (negligible)
    angq < 0 && (angq += 360)
    angq < 135 && (angq += 360)
    return abs(angq - 315)
end

"""
From Sheddy 1968, A General Analytic Solution for Reflection ...
"""
function sharplyboundedR(ea::EigenAngle, M)
    # Unpack
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # Solve equation `D = 0` as a quartic
    q = bookerquartic(ea, M)

    # For high in the ionosphere, we choose 2 solutions that lie close to positive real and
    # negative imaginary axis (315° on the complex plane)
    sort!(q, by=anglefrom315)

    # Constant entries of dispersion matrix `D`
    D12 = M[1,2]
    D32 = M[3,2]
    D33 = C² + M[3,3]

    # Values for two solutions of the Booker quartic corresponding to upgoing waves
    D11_1 = 1 + M[1,1] - q[1]^2
    D13_1 = M[1,3] + q[1]*S
    D31_1 = M[3,1] + q[1]*S

    Δ_1 = D11_1*D33 - D13_1*D31_1
    P_1 = (-D12*D33 + D13_1*D32)/Δ_1
    T_1 = q[1]*P_1

    D11_2 = 1 + M[1,1] - q[2]^2
    D13_2 = M[1,3] + q[2]*S
    D31_2 = M[3,1] + q[2]*S

    Δ_2 = D11_2*D33 - D13_2*D31_2
    P_2 = (-D12*D33 + D13_2*D32)/Δ_2
    T_2 = q[2]*P_2

    # Computation of entries of reflection matrix `R`
    Δ = (T_1*C + P_1)*(C + q[2]) - (T_2*C + P_2)*(C + q[1])
    invΔ = 1/Δ

    R11 = ((T_1*C - P_1)*(C + q[2]) - (T_2*C - P_2)*(C + q[1]))*invΔ  # ∥R∥
    R22 = ((T_1*C + P_1)*(C - q[2]) - (T_2*C + P_2)*(C - q[1]))*invΔ  # ⟂R⟂
    R12 = -2*C*(T_1*P_2 - T_2*P_1)*invΔ  # ⟂R∥
    R21 = -2*C*(q[1] - q[2])*invΔ  # ∥R⟂

    return SMatrix{2,2}(R11, R21, R12, R22)
end

"""
Calculate the derivative of the reflection matrix `R` wrt height `z`.

Expects S to be a tuple of `(S11, S12, S21, S22)`.

The factor of `k` appears explicitly because Budden 1955 derives R′ wrt a height
variable ``s'' which includes k.
"""
function dRdz(dR, R, params, z)
    ω, k = params.ω, params.k
    ea = params.ea
    z0, species, bfield = params.referenceheight, params.species, params.bfield

    M = susceptibility(ω, z0, z, species, bfield)
    S = smatrix(ea, M)

    return -im*k/2*(S[3] + S[4]*R - R*S[1] - R*S[2]*R)
end

function integratethroughionosphere(
    ea::EigenAngle,
    source::AbstractSource,
    fromheight,
    toheight,
    referenceheight,
    species,
    bfield::BField
)
    # Unpack
    ω, k, λ = source.ω, source.k, source.λ  # k in km

    M = susceptibility(ω, referenceheight, fromheight, species, bfield)
    R0 = sharplyboundedR(ea, M)

    params = (ω=ω, k=k, ea=ea, referenceheight=referenceheight, species=species, bfield=bfield)
    prob = ODEProblem(dRdz, R0, (fromheight, toheight), params)
    sol = solve(prob, reltol=1e-6)#, dtmax=λ/20)
end



# Values taken from the end of a random homogeneous exponential LWPC run
θ = complex(65.2520447, -1.5052794)
ω = 150796.4531250  # rad/s
freq = ω/2π

ea = EigenAngle(θ)
source = Source(freq)

Bmag = 0.5172267e-4  # (T) LWPC is actually using Gauss, despite their claim of W/m²=T
dcl = -0.2664399
dcm = -0.2850476
dcn = -0.9207376
bfield = BField(Bmag, dcl, dcm, dcn)

species = Constituent(-fundamentalcharge, mₑ,
                        h -> waitprofile(h, 75, 0.3), collisionprofile)

referenceheight = 50.0
fromheight = 90.3906250
toheight = 44.2968750
height = 70.0
