#==
Functions related to identifying resonant modes ("eigenangles") within the
earth-ionosphere waveguide.

These functions use Budden's model, which represents the ground and ionosphere
as sharply reflecting boundaries. The middle of the waveguide is filled with a
fictitious medium with a refractive index that mimicks the propagation of the
radio wave through free space over curved earth.
==#

# TODO: G as its own type contained within this?
@with_kw struct SharpBoundaryVariables{T<:Number} @deftype T
    G12
    G32
    G33
    G11₁
    G13₁
    G31₁
    Δ₁
    Δ₁⁻¹
    P₁
    T₁
    G11₂
    G13₂
    G31₂
    Δ₂
    Δ₂⁻¹
    P₂
    T₂
    Δ
    Δ⁻¹
end

function Base.show(io::IO, ::MIME"text/plain", s::SharpBoundaryVariables{T}) where {T}
    println(io, "SharpBoundaryVariables{$T}: ")
    for n in fieldnames(SharpBoundaryVariables)
        if n != last(fieldnames(SharpBoundaryVariables))
            println(io, " $n: ", getfield(s, n))
        else
            # To avoid double newlines (one is automatic), use `print` on last field
            print(io, " $n: ", getfield(s, n))
        end
    end
end

@with_kw struct FresnelReflectionVariables{T<:Number} @deftype T
    Ng²
    CNg²
    sqrtNg²mS²
    Rg::SDiagonal{2,T}
end

##########
# Reflection coefficients
##########

function _sharpboundaryreflection(ea::EigenAngle, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    # XXX: `bookerquartic` (really `roots!`) dominates this functions runtime
    bookerquartic!(ea, M)

    # We choose the 2 roots corresponding to upward travelling waves as being
    # those that lie close to the positive real axis and negative imaginary axis
    sortquarticroots!(BOOKER_QUARTIC_ROOTS)

    #==
    Γ = [0 -q 0;
         q 0 -S;
         0 S 0]
    G = Γ² + I + M
    G = [1-q[i]^2+M[1,1] M[1,2] S*q[i]+M[1,3];
         M[2,1] 1-q[i]^2-S²+M[2,2] M[2,3];
         S*q[i]+M[3,1] M[3,2] C²+M[3,3]]
    ==#

    #BUG: Need to check that Sheddy actually wants roots 2, 1 in that order rather
    # than the order used by Pitteway.

    # Precompute
    q = BOOKER_QUARTIC_ROOTS
    q₁S = q[1]*S
    q₂S = q[2]*S

    M11p1 = 1 + M[1,1]

    # Constant entries of dispersion matrix `G`
    G12 = M[1,2]
    G32 = M[3,2]
    G33 = C² + M[3,3]

    G12G33 = G12*G33

    # Values for solutions of the Booker quartic corresponding to upgoing waves
    G11₁ = M11p1 - q[1]^2
    G13₁ = M[1,3] + q₁S
    G31₁ = M[3,1] + q₁S

    Δ₁ = G11₁*G33 - G13₁*G31₁
    Δ₁⁻¹ = 1/Δ₁
    P₁ = (-G12G33 + G13₁*G32)*Δ₁⁻¹
    T₁ = q[1]*P₁ - S*(-G11₁*G32 + G12*G31₁)*Δ₁⁻¹

    G11₂ = M11p1 - q[2]^2
    G13₂ = M[1,3] + q₂S
    G31₂ = M[3,1] + q₂S

    Δ₂ = G11₂*G33 - G13₂*G31₂
    Δ₂⁻¹ = 1/Δ₂
    P₂ = (-G12G33 + G13₂*G32)*Δ₂⁻¹
    T₂ = q[2]*P₂ - S*(-G11₂*G32 + G12*G31₂)*Δ₂⁻¹

    Δ = (T₁*C + P₁)*(C + q[2]) - (T₂*C + P₂)*(C + q[1])
    Δ⁻¹ = 1/Δ

    return SharpBoundaryVariables(G12, G32, G33, G11₁, G13₁, G31₁, Δ₁, Δ₁⁻¹, P₁, T₁,
                                  G11₂, G13₂, G31₂, Δ₂, Δ₂⁻¹, P₂, T₂, Δ, Δ⁻¹)
end

"""
    sharplyboundedreflection(ea, M)

Return ionosphere reflection matrix `R` for a sharply bounded anisotropic ionosphere.

[^Sheddy1968a] introduces the matrix equation ``G E = 0`` in the coordinate system of Budden
where the dispersion matrix ``G = I + M + L``. Nontrivial solutions to the equation require
``det(G) = 0``, which may be written as a quartic in ``q = n cos(θ)``. Elements of the
dispersion matrix are then used to calculate the reflection matrix `R`.

The reflection coefficient matrix for the sharply bounded case is used as a starting solution
for integration of the reflection coefficient matrix through the ionosphere.

# References

[^Sheddy1968a]: C. H. Sheddy, “A General Analytic Solution for Reflection From a Sharply Bounded Anisotropic Ionosphere,” Radio Science, vol. 3, no. 8, pp. 792–795, Aug. 1968.
"""
function sharpboundaryreflection(ea::EigenAngle, M)
    C = ea.cosθ
    C2 = 2*C

    @unpack P₁, T₁, P₂, T₂, Δ⁻¹ = _sharpboundaryreflection(ea, M)

    # NOTE: `BOOKER_QUARTIC_ROOTS` is sorted in _sharpboundaryreflection
    # Could write an `issorted` function, but sorting time is dominated by
    # `upgoing`, which would presumably be required by an `issorted`
    q = BOOKER_QUARTIC_ROOTS

    T₁C = T₁*C
    T₂C = T₂*C

    R11 = ((T₁C - P₁)*(C + q[2]) - (T₂C - P₂)*(C + q[1]))*Δ⁻¹  # ∥R∥
    R22 = ((T₁C + P₁)*(C - q[2]) - (T₂C + P₂)*(C - q[1]))*Δ⁻¹  # ⟂R⟂
    R12 = -C2*(T₁*P₂ - T₂*P₁)*Δ⁻¹  # ⟂R∥
    R21 = -C2*(q[1] - q[2])*Δ⁻¹  # ∥R⟂

    return SMatrix{2,2}(R11, R21, R12, R22)
end

# TODO: Autodiff with Zygote?
function sharpboundaryreflection(ea::EigenAngle, M, ::Derivative_dθ)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    @unpack G12, G32, G33, G11₁, G13₁, G31₁, Δ₁, Δ₁⁻¹, P₁, T₁,
        G11₂, G13₂, G31₂, Δ₂, Δ₂⁻¹, P₂, T₂, Δ, Δ⁻¹ = _sharpboundaryreflection(ea, M)
    q, B = BOOKER_QUARTIC_ROOTS, BOOKER_QUARTIC_COEFFS

    # Precompute some variables
    C2 = 2*C

    T₁C = T₁*C
    T₂C = T₂*C

    M13pM31 = M[1,3] + M[3,1]

    # Additional calculations required for dR/dθ
    dS = C
    dC = -S
    dC² = -S*C2
    dB3 = dS*M13pM31
    dB2 = -dC²*(2 + M[1,1] + M[3,3])
    dB1 = dS/S*B[2] - S*dC²*M13pM31
    dB0 = dC²*(2*C²*(1 + M[1,1]) + M[3,3] + M[2,2] + M[1,1]*(M[3,3] + M[2,2]) -
            M[1,3]*M[3,1] - M[1,2]*M[2,1])

    dq_1 = -(((dB3*q[1] + dB2)*q[1] + dB1)*q[1] + dB0) /
            (((4*B[5]*q[1] + 3*B[4])*q[1] + 2*B[3])*q[1] + B[2])
    dq_2 = -(((dB3*q[2] + dB2)*q[2] + dB1)*q[2] + dB0) /
            (((4*B[5]*q[2] + 3*B[4])*q[2] + 2*B[3])*q[2] + B[2])

    dG33 = dC²

    dG11₁ = -2*q[1]*dq_1
    dG13₁ = dq_1*S + q[1]*dS
    dG31₁ = dG13₁  # dq_1*S + q[1]*dS

    dΔ₁ = dG11₁*G33 + G11₁*dG33 - dG13₁*G31₁ - G13₁*dG31₁
    dΔ₁⁻¹ = -dΔ₁/Δ₁^2

    dP₁ = (-G12*dG33 + dG13₁*G32)*Δ₁⁻¹ + (G13₁*G32 - G12*G33)*dΔ₁⁻¹
    dT₁ = dq_1*P₁ + q[1]*dP₁ -
        dS*(-G11₁*G32 + G12*G31₁)*Δ₁⁻¹ -
        S*(-dG11₁*G32 + G12*dG31₁)*Δ₁⁻¹ -
        S*(-G11₁*G32 + G12*G31₁)*dΔ₁⁻¹

    dG11₂ = -2*q[2]*dq_2
    dG13₂ = dq_2*S + q[2]*dS
    dG31₂ = dG13₂  # dq_2*S + q[2]*dS

    dΔ₂ = dG11₂*G33 + G11₂*dG33 - dG13₂*G31₂ - G13₂*dG31₂
    dΔ₂⁻¹ = -dΔ₂/Δ₂^2

    dP₂ = (-G12*dG33 + dG13₂*G32)*Δ₂⁻¹ + (G13₂*G32 - G12*G33)*dΔ₂⁻¹
    dT₂ = dq_2*P₂ + q[2]*dP₂ -
        dS*(-G11₂*G32 + G12*G31₂)*Δ₂⁻¹ -
        S*(-dG11₂*G32 + G12*dG31₂)*Δ₂⁻¹ -
        S*(-G11₂*G32 + G12*G31₂)*dΔ₂⁻¹

    dΔ = dT₁*C² + T₁*dC² + dT₁*C*q[2] + T₁*dC*q[2] + T₁C*dq_2 + dP₁*C + P₁*dC +
            dP₁*q[2] + P₁*dq_2 - (dT₂*C² + T₂*dC²) -
            (dT₂*C*q[1] + T₂*dC*q[1] + T₂C*dq_1) -
            (dP₂*C + P₂*dC) - (dP₂*q[1] + P₂*dq_1)
    dΔ⁻¹ = -dΔ/Δ^2

    # R
    R11 = ((T₁C - P₁)*(C + q[2]) - (T₂C - P₂)*(C + q[1]))*Δ⁻¹  # ∥R∥
    R22 = ((T₁C + P₁)*(C - q[2]) - (T₂C + P₂)*(C - q[1]))*Δ⁻¹  # ⟂R⟂
    R12 = -C2*(T₁*P₂ - T₂*P₁)*Δ⁻¹  # ⟂R∥
    R21 = -C2*(q[1] - q[2])*Δ⁻¹  # ∥R⟂

    # dR
    dR11 = dΔ⁻¹*R11*Δ +
        Δ⁻¹*(C²*(dT₁ - dT₂) + dC²*(T₁ - T₂) + dC*(T₁*q[2] - P₁ - T₂*q[1] + P₂) +
            C*(dT₁*q[2] + T₁*dq_2 + dP₂ - dP₁ - dT₂*q[1] - T₂*dq_1) +
            dP₂*q[1] + P₂*dq_1 - dP₁*q[2] - P₁*dq_2)
    dR12 = dΔ⁻¹*R12*Δ + dC/C*R12 - C2*(dT₁*P₂ + T₁*dP₂ - dT₂*P₁ - T₂*dP₁)*Δ⁻¹
    dR21 = dΔ⁻¹*R21*Δ + dC/C*R21 - C2*(dq_1 - dq_2)*Δ⁻¹
    dR22 = dΔ⁻¹*R22*Δ +
            Δ⁻¹*(C²*(dT₁ - dT₂) + dC²*(T₁ - T₂) + dC*(T₂*q[1] - T₁*q[2] + P₁ - P₂) +
            C*(dP₁ - dT₁*q[2] + dT₂*q[1] - T₁*dq_2 + T₂*dq_1 - dP₂) +
            dP₂*q[1] + P₂*dq_1 - dP₁*q[2] - P₁*dq_2)

    return @SMatrix [R11 R12;
                     R21 R22;
                     dR11 dR12;
                     dR21 dR22]
end

"""
    wmatrix(ea, T)

Compute submatrix elements of `W` for solving the differential equation of reflection matrix
`R` wrt `z`.

``W`` is also known as ``S`` in many texts.

Following Budden's [^Budden1955a] formalism for the reflection matrix of a plane wave
obliquely incident on the ionosphere, the wave below the ionosphere can be resolved into
upgoing and downgoing waves of elliptical polarization, each of whose components are
themselves resolved into a component with the electric field in the plane of propagation and
a component perpendicular to the plane of propagation. The total field can be written in
matrix form as ``e = L q`` where ``L`` is a 4×4 matrix that simply selects and specifies the
incident angle of the components and ``q`` is a column matrix of the complex amplitudes of
the component waves. By inversion, ``q = L⁻¹ e`` and its derivative wrt height `z` is
``q′ = -i L⁻¹ T L q = -½ i S q``. Then ``S = 2 L⁻¹ T L``, which is calculated by this
function, describes the change in amplitude of the upgoing and downgoing component waves.

See also: [`tmatrix`](@ref), [`susceptibility`](@ref)

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
[^Sheddy1968]: C. H. Sheddy, R. Pappert, Y. Gough, and W. Moler, “A Fortran program for mode constants in an earth-ionosphere waveguide,” Naval Electronics Laboratory Center, San Diego, CA, Interim Report 683, May 1968.
"""
function wmatrix(ea::EigenAngle, T::TMatrix)
    C = ea.cosθ
    Cinv = ea.secθ

    # Precompute
    T12Cinv = T[1,2]*Cinv
    T14Cinv = T[1,4]*Cinv
    T32Cinv = T[3,2]*Cinv
    T34Cinv = T[3,4]*Cinv
    CT41 = C*T[4,1]

    #==
    Sheddy 1968 has at least 1 "typo": S11b should be written T14/C+CT41
    Budden 1988 has at least 1 typo: -T11 - T14 + CT41 + T44 for element [1,1] of [2,1]

    The analytical solution below was checked using SymPy.
    There is also a test against 2*inv(L)*T*L:

    T = SMatrix{4,4}(T11, 0, T31, T41,
                     T12, 0, T32, T42,
                     0, 1, 0, 0,
                     T14, 0, T34, T44)
    L = SMatrix{4,4}(C, 0, 0, 1,
                     0, -1, -C, 0,
                     -C, 0, 0, 1,
                     0, -1, C, 0)
    # 2*inv(L) = Linv2
    L2inv = SMatrix{4,4}(Cinv, 0, -Cinv, 0,
                        0, -1, 0, -1,
                        0, -Cinv, 0, Cinv,
                        1, 0, 1, 0)
    W = Linv2*T*L

    # W = 2*(L\T)*L

    ---

    W = | W11 | W12 |
        | W21 | W22 |

    W11 = | a11+a11r | -b11 |
          | -c11     | d11  |

    W12 = | a21+a21r | -b11 |
          | c12      | d12  |

    W21 = | a21-a21r | b21  |
          | c11      | -d12 |

    W22 = | a11-a11r | b21  |
          | -c12     | -d11 |
    ==#

    a11 = T[1,1] + T[4,4]
    a11r = T14Cinv + CT41
    a21 = T[4,4] - T[1,1]
    a21r = T14Cinv - CT41

    b11 = T12Cinv + T[4,2]
    b21 = T12Cinv - T[4,2]

    c11 = T[3,1] + T34Cinv
    c12 = T[3,1] - T34Cinv

    d11 = C + T32Cinv
    d12 = T32Cinv - C

    # Form the four 2x2 submatrices of `S`
    # Ttype is already promoted type of ea and M
    W11 = SMatrix{2,2,ComplexF64,4}(a11+a11r, -c11, -b11, d11)
    W12 = SMatrix{2,2,ComplexF64,4}(a21+a21r, c12, -b11, d12)
    W21 = SMatrix{2,2,ComplexF64,4}(a21-a21r, c11, b21, -d12)
    W22 = SMatrix{2,2,ComplexF64,4}(a11-a11r, -c12, b21, -d11)

    return W11, W21, W12, W22
end

"""
    dwmatrixdθ(ea, M, T)

Return the 4 submatrix elements of the derivative of the `W` matrix wrt θ.

See also: [`wmatrix`](@ref), [`susceptibility`](@ref), [`tmatrix`](@ref)
"""
function dwmatrixdθ(ea::EigenAngle, M, T::TMatrix)
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = ea.secθ
    C²inv = Cinv^2

    dS = C
    dC = -S
    dC² = -2*S*C
    dCinv = S*C²inv

    den = inv(1 + M[3,3])

    dT11 = -dS*M[3,1]*den
    dT12 = dS*M[3,2]*den
    # dT13 = 0
    dT14 = dC²*den
    # dT21 = 0
    # dT22 = 0
    # dT23 = 0
    # dT24 = 0
    # dT31 = 0
    dT32 = dC²
    # dT33 = 0
    dT34 = dS*M[2,3]*den
    # dT41 = 0
    # dT42 = 0
    # dT43 = 0
    dT44 = -dS*M[1,3]*den

    dt12Cinv = dT12*Cinv + T[1,2]*dCinv
    dt14Cinv = dT14*Cinv + T[1,4]*dCinv
    dt32Cinv = dT32*Cinv + T[3,2]*dCinv
    dt34Cinv = dT34*Cinv + T[3,4]*dCinv
    dt41C = dC*T[4,1]

    ds11a = dT11 + dT44
    dd11a = dT11 - dT44
    ds11b = dt14Cinv + dt41C
    dd11b = dt14Cinv - dt41C
    ds12 = dt12Cinv
    dd12 = dt12Cinv
    ds21 = dt34Cinv
    dd21 = -dt34Cinv
    ds22 = dC + dt32Cinv
    dd22 = dC - dt32Cinv

    # Form the four 2x2 submatrices of `dS`
    dW11 = SMatrix{2,2,ComplexF64,4}(ds11a+ds11b, -ds21, -ds12, ds22)
    dW12 = SMatrix{2,2,ComplexF64,4}(-dd11a+dd11b, dd21, -ds12, -dd22)
    dW21 = SMatrix{2,2,ComplexF64,4}(-dd11a-dd11b, ds21, dd12, dd22)
    dW22 = SMatrix{2,2,ComplexF64,4}(ds11a-ds11b, -dd21, dd12, -ds22)

    return dW11, dW21, dW12, dW22
end

"""
    dRdz(R, params, z)

Return the differential of the reflection matrix `R` wrt height `z`.

Following the Budden [^Budden1955a] formalism for the reflection of an (obliquely) incident
plane wave from a horizontally stratified ionosphere, the differential of the reflection
matrix `R` with height `z` can be described by ``2i R′ = W⁽²¹⁾ + W⁽²²⁾R - RW⁽¹¹⁾ - RW⁽¹²⁾R``.
Integrating ``R′`` wrt height `z`, gives the reflection matrix ``R`` for the ionosphere.

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955.
"""
function dRdz(R, params, z)
    ea, frequency, waveguide = params
    @unpack bfield, species = waveguide

    k = frequency.k

    # XXX: `susceptibility` takes as much time as `tmatrix` and `wmatrix` combined
    # TODO: Maybe use a function approximator version of M b/c it's not dependent on θ
    M = susceptibility(z, frequency, bfield, species)
    T = tmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)

    # TEMP BUG TODO XXX
    # return R*1.01

    # the factor k/(2i) isn't explicitly in Budden1955a b/c of his change of variable ``s = kz``
    return -im/2*k*(W21 + W22*R - R*W11 - R*W12*R)
end

function dRdθdz(RdRdθ, params, z)
    ea, frequency, waveguide = params
    @unpack bfield, species = waveguide

    k = frequency.k

    M = susceptibility(z, frequency, bfield, species)
    T = tmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)
    dW11, dW21, dW12, dW22 = dwmatrixdθ(ea, M, T)

    R = RdRdθ[SVector(1,2),:]
    dRdθ = RdRdθ[SVector(3,4),:]

    dz = -im/2*k*(W21 + W22*R - R*W11 - R*W12*R)
    dθdz = -im/2*k*(dW21 + dW22*R + W22*dRdθ - (dRdθ*W11 + R*dW11) -
                    (dRdθ*W12*R + R*dW12*R + R*W12*dRdθ))

    return vcat(dz, dθdz)
end

function integratedreflection(ea::EigenAngle, frequency::Frequency, waveguide::HomogeneousWaveguide{T}) where T
    @unpack bfield, species = waveguide

    Mtop = susceptibility(TOPHEIGHT, frequency, bfield, species)
    Rtop = sharpboundaryreflection(ea, Mtop)

    # TODO: Pass parameters with tolerances, integration method
    # TODO: Does saving actually alter performance?
    prob = ODEProblem{false}(dRdz, Rtop, (TOPHEIGHT, BOTTOMHEIGHT), (ea, frequency, waveguide))

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8,
                save_on=false, save_start=false, save_end=true)

    R::typeof(Rtop) = sol[end]  # TODO: why is sol::Any?

    return R
end

# This is kept as a completely separate function because the size of the matrix being
# integrated is different and therefore the size of sol[end] is different too
# The derivative terms are intertwined with the non-derivative terms so we can't do only
# the derivative terms
function integratedreflection(ea::EigenAngle, frequency::Frequency, waveguide::HomogeneousWaveguide{T}, ::Derivative_dθ) where T
    @unpack bfield, species = waveguide

    Mtop = susceptibility(TOPHEIGHT, frequency, bfield, species)

    RdRdθtop = sharpboundaryreflection(ea, Mtop, Derivative_dθ())

    prob = ODEProblem{false}(dRdθdz, RdRdθtop, (TOPHEIGHT, BOTTOMHEIGHT), (ea, frequency, waveguide))

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8,
                save_on=false, save_start=false, save_end=true)

    R::typeof(RdRdθtop) = sol[end]  # TODO: why is sol::Any?

    return R
end


##########
# Ground reflection coefficient matrix
##########

function _fresnelreflection(ea, ground, frequency)
    C, S² = ea.cosθ, ea.sin²θ
    ω = frequency.ω

    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*E0))

    CNg² = C*Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²mS²)/(CNg² + sqrtNg²mS²)
    Rg22 = (C - sqrtNg²mS²)/(C + sqrtNg²mS²)

    Rg = SDiagonal(Rg11, Rg22)

    return FresnelReflectionVariables(Ng², CNg², sqrtNg²mS², Rg)
end

"""
    fresnelreflection(ea, ground, frequency)

Return the Fresnel reflection coefficient matrix for the ground free-space interface at the
ground (``z = 0``). Follows the formulation in [^Morfitt1976] pages 25-26.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency)
    @unpack Rg = _fresnelreflection(ea, ground, frequency)

    return Rg
end

function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency, ::Derivative_dθ)
    C, S, S² = ea.cosθ, ea.sinθ, ea.sin²θ
    ω = frequency.ω

    @unpack Rg, Ng², CNg², sqrtNg²mS² = _fresnelreflection(ea, ground, frequency)

    S2 = 2*S

    dRg11 = (S2*Ng²*(1 - Ng²))/(sqrtNg²mS²*(CNg² + sqrtNg²mS²)^2)
    dRg22 = (S2*(C - sqrtNg²mS²))/(sqrtNg²mS²*(sqrtNg²mS² + C))

    dRg = SDiagonal(dRg11, dRg22)

    return Rg, dRg
end

##########
# Identify EigenAngles
##########

"""
    modalequation(R, Rg)

Return value of the determinental mode equation ``det(Rg R - I)`` given `R` and `Rg`.

If the reflection matrix ``R₁`` represents the reflection coefficient at height ``z₁``, then
the reflection coefficient at ``z₂`` is ``R₂ = R₁ exp(2ik(z₂ - z₁) sin(θ))``. For a wave that
starts at the ground, reflects off the ionosphere at level ``h`` with ``R``, and then at the
ground with ``Rg``, the fields must be multiplied by ``̅Rg R exp(-2ikh sin(θ))``. For a
self-consistent waveguide mode, the resulting wave must be identifcal with the original
wave, and a necessary and sufficient condition for this is that one eigenvalue of this matrix
is unity. This leads to the mode equation [^Budden1962].

# References

[^Budden1962]: K. G. Budden and N. F. Mott, “The influence of the earth’s magnetic field on
radio propagation by wave-guide modes,” Proceedings of the Royal Society of London. Series A.
Mathematical and Physical Sciences, vol. 265, no. 1323, pp. 538–553, Feb. 1962.
"""
function modalequation(R, Rg)
    return det(Rg*R - I)
end

"""
See https://folk.ntnu.no/hanche/notes/diffdet/diffdet.pdf
"""
function modalequationdθ(R, dR, Rg, dRg)
    A = Rg*R - I
    dA = dRg*R + Rg*dR
    return det(A)*tr(A\dA)  # much faster than A\dA TODO: preferred approach?
end

function solvemodalequation(ea::EigenAngle, frequency::Frequency, waveguide::HomogeneousWaveguide{T}) where T
    R = integratedreflection(ea, frequency, waveguide)
    Rg = fresnelreflection(ea, waveguide.ground, frequency)

    f = modalequation(R, Rg)
    return f
end

"""
This returns R and Rg in addition to df because the only time this function is needed, we also
need R and Rg (in excitationfactors).
"""
function solvemodalequationdθ(ea::EigenAngle, frequency::Frequency, waveguide::HomogeneousWaveguide{T}) where T
    RdR = integratedreflection(ea, frequency, waveguide, Derivative_dθ())
    R = RdR[SVector(1,2),:]
    dR = RdR[SVector(3,4),:]

    Rg, dRg = fresnelreflection(ea, waveguide.ground, frequency, Derivative_dθ())

    df = modalequationdθ(R, dR, Rg, dRg)
    return df, R, Rg
end

"""
`tolerance` is how close solution of modal equation gets to 0. Difference in θ between
`tolerance=1e-6` and `tolerance=1e-8` is less than 1e-5° in both real and imaginary
components.
"""
function findmodes(origcoords::AbstractVector, frequency::Frequency, waveguide::HomogeneousWaveguide{T}, tolerance=1e-8) where {T}

    # WARNING: If tolerance of mode finder is much less than the R integration
    # tolerance, it may possible multiple identical modes will be identified.

    est_num_nodes = ceil(Int, length(origcoords)*1.1)

    zroots, zpoles = grpf(z->solvemodalequation(EigenAngle(z), frequency, waveguide),
                          origcoords, GRPFParams(est_num_nodes, tolerance, true))

    return EigenAngle.(zroots)
end
