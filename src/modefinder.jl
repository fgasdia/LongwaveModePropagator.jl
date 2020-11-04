#==
Functions related to identifying resonant modes ("eigenangles") within the
earth-ionosphere waveguide.

These functions use Budden's model, which represents the ground and ionosphere
as sharply reflecting boundaries. The middle of the waveguide is filled with a
fictitious medium with a refractive index that mimicks the propagation of the
radio wave through free space over curved earth.
==#

abstract type ModeEquation end

struct PhysicalModeEquation{W<:HomogeneousWaveguide} <: ModeEquation
    ea::EigenAngle
    frequency::Frequency
    waveguide::W
end
PhysicalModeEquation(f::Frequency, w::HomogeneousWaveguide) =
    PhysicalModeEquation(EigenAngle(complex(0.0)), f, w)
setea(ea::EigenAngle, p::PhysicalModeEquation) = PhysicalModeEquation(ea, p.frequency, p.waveguide)

struct IntegrationParams{T}
    tolerance::Float64
    solver::T
    force_dtmin::Bool
end
const DEFAULT_INTEGRATIONPARAMS = IntegrationParams(1e-6, BS5(), false)

const GRPF_PARAMS = Ref{GRPFParams}()
get_grpf_params() = GRPF_PARAMS[]
set_grpf_params() = GRPF_PARAMS[] = GRPFParams(100000, 1e-5, true)
set_grpf_params(p::GRPFParams) = GRPF_PARAMS[] = p


"""
    isroot(z::Complex; atol=1e-2)

Return `true` if both real and imaginary components of `z` are approximately equal to 0.
"""
function isroot(z::Complex; atol=1e-2)
    rz, iz = reim(z)
    isapprox(rz, 0, atol=atol) && isapprox(iz, 0, atol=atol)
end

"""
    isroot(x::Real; atol=1e-2)

Return `true` if the value of `x` is approximately equal to 0.
"""
isroot(x::Real; atol=1e-2) = isapprox(z, 0, atol=atol)

"""
    filterroots!(roots, frequency, waveguide; atol=1e-2)

Remove elements from `roots` if they are not valid roots (determined by `isroot`) to the modal
equation.

The absolute tolerance `atol` is passed through to `isroot`.
"""
function filterroots!(roots, frequency, waveguide; atol=1e-2)
    modeequation = PhysicalModeEquation(EigenAngle(0.0), frequency, waveguide)
    filterroots!(roots, modeequation, atol=atol)
end

"""
    filterroots!(roots, modeequation::PhysicalModeEquation; atol=1e-2)

Use existing `modeequation` to validate `roots`.

!!! note

    `modeequation.ea` will be modified by this function.
"""
function filterroots!(roots, modeequation::PhysicalModeEquation; atol=1e-2)
    i = 1
    while i <= length(roots)
        modeequation = setea(EigenAngle(roots[i]), modeequation)
        f = solvemodalequation(modeequation)
        isroot(f, atol=atol) ? (i += 1) : deleteat!(roots, i)
    end
    return nothing
end


##########
# Reflection coefficients
##########

"""
    sharpboundaryreflection(ea, M)

Return ionosphere reflection matrix `R` for a sharply bounded anisotropic ionosphere
referenced to the ground.

The reflection coefficient matrix for the sharply bounded case is used as a starting solution
for integration of the reflection coefficient matrix through the ionosphere.

From the differential equations for fields in the ionosphere ``de/dz = -ikTe``, assuming
a homogeneous ionosphere ``(T - q)e = 0``. This only has a non-trivial solution if
``det(T - q) = 0``, which is a fourth degree equation in ``q``, the Booker quartic.

From linear algebra theory, the four ``q``s are eigenvalues of ``T`` and the ``e``s are the
associated eigencolumns, then ``Teᵢ = qᵢeᵢ``. If any ``eᵢ`` is multiplied by a constant, it
remains an eigencolumn of ``T``.

The reflection coefficient is the ratio of field components, thus we are free to choose any
consistent scaling of the fields ``e``. This function uses the scaling in [^Budden1988]
pg. 190 and the expression for the reflection coefficient in [^Budden1988] pg. 307.

# References

[^Sheddy1968a]: C. H. Sheddy, “A General Analytic Solution for Reflection From a Sharply
Bounded Anisotropic Ionosphere,” Radio Science, vol. 3, no. 8, pp. 792–795, Aug. 1968.
[^Budden1988]:
"""
function sharpboundaryreflection(ea::EigenAngle, M)
    T = tmatrix(ea, M)
    q, B = bookerquartic(T)

    sort!(q, by=upgoing)

    C = ea.cosθ

    # Compute fields using the ratios from [^Budden1988] pg 190
    a1 = -(T[1,1] + T[4,4])
    a2 = T[1,1]*T[4,4] - T[1,4]*T[4,1]
    a3 = T[1,2]
    a4 = T[1,4]*T[4,2] - T[1,2]*T[4,4]
    a5 = T[4,2]
    a6 = T[1,2]*T[4,1] - T[1,1]*T[4,2]

    e = MMatrix{4,2,eltype(T),8}(undef)
    @inbounds for i = 1:2
        A = q[i]^2 + a1*q[i] + a2

        e[1,i] = a3*q[i] + a4
        e[2,i] = A
        e[3,i] = q[i]*A
        e[4,i] = a5*q[i] + a6
    end

    # Compute reflection matrix referenced to z=0 [^Budden1988] pg 307
    d = SMatrix{2,2}(C*e[4,1]-e[1,1], -C*e[2,1]+e[3,1], C*e[4,2]-e[1,2], -C*e[2,2]+e[3,2])
    u = SMatrix{2,2}(C*e[4,1]+e[1,1], -C*e[2,1]-e[3,1], C*e[4,2]+e[1,2], -C*e[2,2]-e[3,2])

    R = d/u

    return R
end

function sharpboundaryreflection(ea::EigenAngle, M, ::Dθ)
    S, C = ea.sinθ, ea.cosθ

    T = tmatrix(ea, M)
    dT = tmatrix(ea, M, Dθ())
    q, B = bookerquartic(ea, M)
    sort!(q, by=upgoing)
    dq = bookerquartic(ea, M, q, B, Dθ())

    # Compute fields using the ratios from [^Budden1988] pg 190
    a1 = -(T[1,1] + T[4,4])
    a2 = T[1,1]*T[4,4] - T[1,4]*T[4,1]
    a3 = T[1,2]
    a4 = T[1,4]*T[4,2] - T[1,2]*T[4,4]
    a5 = T[4,2]
    a6 = T[1,2]*T[4,1] - T[1,1]*T[4,2]

    # Many `dT` components are 0
    da1 = -dT[1,1] - dT[4,4]
    da2 = T[1,1]*dT[4,4] - T[4,1]*dT[1,4] + T[4,4]*dT[1,1]
    da3 = dT[1,2]
    da4 = -T[1,2]*dT[4,4] + T[4,2]*dT[1,4] - T[4,4]*dT[1,2]
    # da5 = 0
    da6 = T[4,1]*dT[1,2] - T[4,2]*dT[1,1]

    e = MMatrix{4,2,eltype(T),8}(undef)
    de = MMatrix{4,2,eltype(T),8}(undef)
    @inbounds for i = 1:2
        A = q[i]^2 + a1*q[i] + a2
        dA = a1*dq[i] + q[i]*da1 + 2*q[i]*dq[i] + da2

        e[1,i] = a3*q[i] + a4
        e[2,i] = A
        e[3,i] = q[i]*A
        e[4,i] = a5*q[i] + a6

        de[1,i] = a3*dq[i] + q[i]*da3 + da4
        de[2,i] = dA
        de[3,i] = A*dq[i] + q[i]*dA
        de[4,i] = a5*dq[i] + da6  # + q[i]*da5 = 0
    end

    d = SMatrix{2,2}(C*e[4,1]-e[1,1], -C*e[2,1]+e[3,1], C*e[4,2]-e[1,2], -C*e[2,2]+e[3,2])
    dd = SMatrix{2,2}(-S*e[4,1]+C*de[4,1]-de[1,1], S*e[2,1]-C*de[2,1]+de[3,1],
                      -S*e[4,2]+C*de[4,2]-de[1,2], S*e[2,2]-C*de[2,2]+de[3,2])
    u = SMatrix{2,2}(C*e[4,1]+e[1,1], -C*e[2,1]-e[3,1], C*e[4,2]+e[1,2], -C*e[2,2]-e[3,2])
    du = SMatrix{2,2}(-S*e[4,1]+C*de[4,1]+de[1,1], S*e[2,1]-C*de[2,1]-de[3,1],
                      -S*e[4,2]+C*de[4,2]+de[1,2], S*e[2,2]-C*de[2,2]-de[3,2])

    R = d/u

    den = inv((u[1,1]*u[2,2] - u[1,2]*u[2,1])^2)
    dR11 = (-d[1,1]*u[2,2] + d[1,2]*u[2,1])*(u[1,1]*du[2,2] - u[1,2]*du[2,1] -
        u[2,1]*du[1,2] + u[2,2]*du[1,1]) + (u[1,1]*u[2,2] - u[1,2]*u[2,1])*(d[1,1]*du[2,2] -
        d[1,2]*du[2,1] - u[2,1]*dd[1,2] + u[2,2]*dd[1,1])
    dR21 = (-d[2,1]*u[2,2] + d[2,2]*u[2,1])*(u[1,1]*du[2,2] - u[1,2]*du[2,1] -
        u[2,1]*du[1,2] + u[2,2]*du[1,1]) + (u[1,1]*u[2,2] - u[1,2]*u[2,1])*(d[2,1]*du[2,2] -
        d[2,2]*du[2,1] - u[2,1]*dd[2,2] + u[2,2]*dd[2,1])
    dR12 = (d[1,1]*u[1,2] - d[1,2]*u[1,1])*(u[1,1]*du[2,2] - u[1,2]*du[2,1] -
        u[2,1]*du[1,2] + u[2,2]*du[1,1]) + (u[1,1]*u[2,2] - u[1,2]*u[2,1])*(-d[1,1]*du[1,2] +
        d[1,2]*du[1,1] + u[1,1]*dd[1,2] - u[1,2]*dd[1,1])
    dR22 = (d[2,1]*u[1,2] - d[2,2]*u[1,1])*(u[1,1]*du[2,2] - u[1,2]*du[2,1] -
        u[2,1]*du[1,2] + u[2,2]*du[1,1]) + (u[1,1]*u[2,2] - u[1,2]*u[2,1])*(-d[2,1]*du[1,2] +
        d[2,2]*du[1,1] + u[1,1]*dd[2,2] - u[1,2]*dd[2,1])
    dR = SMatrix{2,2}(dR11*den, dR21*den, dR12*den, dR22*den)

    #==
    [R11 R12;
     R21 R22;
     dR11 dR12;
     dR21 dR22]
    ==#
    return vcat(R, dR)
end

"""
    wmatrix(ea::EigenAngle, T)

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
function wmatrix(ea::EigenAngle, T)
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
    W11 = SMatrix{2,2}(a11+a11r, -c11, -b11, d11)
    W12 = SMatrix{2,2}(a21+a21r, c12, -b11, d12)
    W21 = SMatrix{2,2}(a21-a21r, c11, b21, -d12)
    W22 = SMatrix{2,2}(a11-a11r, -c12, b21, -d11)

    return W11, W21, W12, W22
end

"""
    dwmatrixdθ(ea, M, T)

Return the 4 submatrix elements of the derivative of the `W` matrix wrt θ.

See also: [`wmatrix`](@ref), [`susceptibility`](@ref), [`tmatrix`](@ref)
"""
function dwmatrixdθ(ea::EigenAngle, M, T)
    C, S, C² = ea.cosθ, ea.sinθ, ea.cos²θ
    Cinv = ea.secθ
    C²inv = Cinv^2

    dC = -S
    dCinv = S*C²inv

    dT = tmatrix(ea, M, Dθ())  # many zeros

    dt12Cinv = dT[1,2]*Cinv + T[1,2]*dCinv
    dt14Cinv = dT[1,4]*Cinv + T[1,4]*dCinv
    dt32Cinv = dT[3,2]*Cinv + T[3,2]*dCinv
    dt34Cinv = dT[3,4]*Cinv + T[3,4]*dCinv
    dt41C = dC*T[4,1]

    ds11a = dT[1,1] + dT[4,4]
    dd11a = dT[1,1] - dT[4,4]
    ds11b = dt14Cinv + dt41C
    dd11b = dt14Cinv - dt41C
    ds12 = dt12Cinv
    dd12 = dt12Cinv
    ds21 = dt34Cinv
    dd21 = -dt34Cinv
    ds22 = dC + dt32Cinv
    dd22 = dC - dt32Cinv

    # Form the four 2x2 submatrices of `dW`
    dW11 = SMatrix{2,2}(ds11a+ds11b, -ds21, -ds12, ds22)
    dW12 = SMatrix{2,2}(-dd11a+dd11b, dd21, -ds12, -dd22)
    dW21 = SMatrix{2,2}(-dd11a-dd11b, ds21, dd12, dd22)
    dW22 = SMatrix{2,2}(ds11a-ds11b, -dd21, dd12, -ds22)

    return dW11, dW21, dW12, dW22
end

"""
    dRdz(R, modeequation, z)

Return the differential of the reflection matrix `R` wrt height `z`.

Following the Budden [^Budden1955a] formalism for the reflection of an (obliquely) incident
plane wave from a horizontally stratified ionosphere, the differential of the reflection
matrix `R` with height `z` can be described by ``2i R′ = W⁽²¹⁾ + W⁽²²⁾R - RW⁽¹¹⁾ - RW⁽¹²⁾R``.
Integrating ``R′`` wrt height `z`, gives the reflection matrix ``R`` for the ionosphere.

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing
reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227,
no. 1171, pp. 516–537, Feb. 1955.
"""
function dRdz(R, modeequation, z)
    @unpack ea, frequency, waveguide = modeequation

    k = frequency.k

    M = susceptibility(z, frequency, waveguide)
    T = tmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)

    # the factor k/(2i) isn't explicitly in Budden1955a b/c of his change of variable ``s = kz``
    return -1im/2*k*(W21 + W22*R - R*W11 - R*W12*R)
end

function dRdθdz(RdRdθ, modeequation, z)
    @unpack ea, frequency, waveguide = modeequation

    k = frequency.k

    M = susceptibility(z, frequency, waveguide)
    T = tmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)
    dW11, dW21, dW12, dW22 = dwmatrixdθ(ea, M, T)

    R = RdRdθ[SVector(1,2),:]
    dRdθ = RdRdθ[SVector(3,4),:]

    dz = -1im/2*k*(W21 + W22*R - R*W11 - R*W12*R)
    dθdz = -1im/2*k*(dW21 + dW22*R + W22*dRdθ - (dRdθ*W11 + R*dW11) -
                    (dRdθ*W12*R + R*dW12*R + R*W12*dRdθ))

    return vcat(dz, dθdz)
end

function integratedreflection(modeequation::PhysicalModeEquation,
    params::IntegrationParams=DEFAULT_INTEGRATIONPARAMS) where T

    @unpack tolerance, solver, force_dtmin = params

    Mtop = susceptibility(TOPHEIGHT, modeequation)
    Rtop = sharpboundaryreflection(modeequation.ea, Mtop)

    prob = ODEProblem{false}(dRdz, Rtop, (TOPHEIGHT, BOTTOMHEIGHT), modeequation)

    # timed with test_modefinder.jl parallel_integration()
    #==
    |  tol  |       method      | time, thrds ms | memory MB |           notes        |
    |-------|-------------------|----------------|-----------|------------------------|
    | 1e-8  | Vern8             |  755, 332      |   12.5    | threads dtmin warnings |
    | 1e-6  | Vern8             |                |           |                        |
    | 1e-8  | Tsit5             |  1199, 567     |    8      |                        |
    | 1e-8  | AutoVern7(Rodas5) |  1011, 459     |    37     | threads dtmin warnings |
    | 1e-8  | Vern6             |  972, 508      |    9      |                        |
    | 1e-8  | BS5               |  892, 442      |    9      |                        |
    ==#

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, solver, abstol=tolerance, reltol=tolerance,
                force_dtmin=force_dtmin,
                save_on=false, save_start=false, save_end=true)

    R = sol[end]

    return R
end

# This is kept as a completely separate function because the size of the matrix being
# integrated is different and therefore the size of sol[end] is different too
# The derivative terms are intertwined with the non-derivative terms so we can't do only
# the derivative terms
function integratedreflection(modeequation::PhysicalModeEquation, ::Dθ,
    params::IntegrationParams=DEFAULT_INTEGRATIONPARAMS) where T

    @unpack tolerance, solver, force_dtmin = params

    Mtop = susceptibility(TOPHEIGHT, modeequation)
    RdRdθtop = sharpboundaryreflection(modeequation.ea, Mtop, Dθ())

    prob = ODEProblem{false}(dRdθdz, RdRdθtop, (TOPHEIGHT, BOTTOMHEIGHT), modeequation)

    # NOTE: When save_on=false, don't try interpolating the solution!
    # Purposefully higher tolerance than non-derivative version
    sol = solve(prob, solver, abstol=tolerance, reltol=tolerance,
                force_dtmin=force_dtmin,
                save_on=false, save_start=false, save_end=true)

    R = sol[end]

    return R
end


##########
# Ground reflection coefficient matrix
##########

"""
    fresnelreflection(ea, ground, frequency)

Return the Fresnel reflection coefficient matrix for the ground free-space interface at the
ground (``z = 0``). Follows the formulation in [^Morfitt1976] pages 25-26.

# References

[^Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976.
"""
function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency)
    C, S² = ea.cosθ, ea.sin²θ
    ω = frequency.ω

    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*E0))

    CNg² = C*Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²mS²)/(CNg² + sqrtNg²mS²)
    Rg22 = (C - sqrtNg²mS²)/(C + sqrtNg²mS²)

    Rg = SDiagonal(Rg11, Rg22)

    return Rg
end
fresnelreflection(m::PhysicalModeEquation) =
    fresnelreflection(m.ea, m.waveguide.ground, m.frequency)

function fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency, ::Dθ)
    C, S, S² = ea.cosθ, ea.sinθ, ea.sin²θ
    S2 = 2*S
    ω = frequency.ω

    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*E0))

    CNg² = C*Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²mS²)/(CNg² + sqrtNg²mS²)
    Rg22 = (C - sqrtNg²mS²)/(C + sqrtNg²mS²)

    Rg = SDiagonal(Rg11, Rg22)

    dRg11 = (S2*Ng²*(1 - Ng²))/(sqrtNg²mS²*(CNg² + sqrtNg²mS²)^2)
    dRg22 = (S2*(C - sqrtNg²mS²))/(sqrtNg²mS²*(sqrtNg²mS² + C))

    dRg = SDiagonal(dRg11, dRg22)

    return Rg, dRg
end
fresnelreflection(m::PhysicalModeEquation, ::Dθ) =
    fresnelreflection(m.ea, m.waveguide.ground, m.frequency, Dθ())

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
    return det(A)*tr(A\dA)
end

function solvemodalequation(modeequation::PhysicalModeEquation,
    integrationparams::IntegrationParams=DEFAULT_INTEGRATIONPARAMS)

    R = integratedreflection(modeequation, integrationparams)
    Rg = fresnelreflection(modeequation)

    f = modalequation(R, Rg)
    return f
end
function solvemodalequation(θ, modeequation::PhysicalModeEquation,
    integrationparams::IntegrationParams=DEFAULT_INTEGRATIONPARAMS)
    # Convenience function for `grpf`
    modeequation = setea(EigenAngle(θ), modeequation)
    solvemodalequation(modeequation, integrationparams)
end

"""
This returns R and Rg in addition to df because the only time this function is needed, we also
need R and Rg (in excitationfactors).
"""
function solvemodalequation(modeequation::PhysicalModeEquation, ::Dθ,
    integrationparams::IntegrationParams=DEFAULT_INTEGRATIONPARAMS)

    RdR = integratedreflection(modeequation, Dθ(), integrationparams)
    R = RdR[SVector(1,2),:]
    dR = RdR[SVector(3,4),:]

    Rg, dRg = fresnelreflection(modeequation, Dθ())

    df = modalequationdθ(R, dR, Rg, dRg)
    return df, R, Rg
end

function solvemodalequation(θ, modeequation::PhysicalModeEquation, ::Dθ,
    integrationparams::IntegrationParams=DEFAULT_INTEGRATIONPARAMS)
    # Convenience function for `grpf`
    modeequation = setea(EigenAngle(θ), modeequation)
    solvemodalequation(modeequation, Dθ(), integrationparams)
end

function findmodes(modeequation::ModeEquation, origcoords;
    grpfparams::GRPFParams=get_grpf_params(),
    integrationparams::IntegrationParams=DEFAULT_INTEGRATIONPARAMS)

    # WARNING: If tolerance of mode finder is much less than the R integration
    # tolerance, it may possible multiple identical modes will be identified. Checks for
    # valid and redundant modes help ensure valid eigenangles are returned from this function.

    roots, poles = grpf(θ->solvemodalequation(θ, modeequation, integrationparams),
                        origcoords, grpfparams)

    # Ensure roots are valid solutions to the modal equation
    filterroots!(roots, modeequation)

    # Remove any redundant modes
    # if tolerance is 1e-8, this rounds to 7 decimal places
    sort!(roots, by=reim)
    ndigits = round(Int, abs(log10(grpfparams.tolerance)+1), RoundDown)
    unique!(z->round(z, digits=ndigits), roots)

    return EigenAngle.(roots)
end

"""
    triangulardomain(Za, Zb, Zc, Δr)

Generate initial mesh node coordinates for a grid-aligned right triangle domain
∈ {`Za`, `Zb`, `Zc`} with initial mesh step `Δr`.

This function generates a mesh grid for particular right triangles where two
sides of the triangle are aligned to the underlying real/imaginary grid. Examples
of such triangles are:

a ----- b
|     /
|   /
| /
c

where
    - a, b have greatest extent in x
    - a, c have greatest extent in y
"""
function triangulardomain(Za::Complex, Zb::Complex, Zc::Complex, Δr)
    rZa, iZa = reim(Za)
    rZb, iZb = reim(Zb)
    rZc, iZc = reim(Zc)

    #==
    # Check if this is a right triangle
    validtriangle = true
    if rZa == rZb == rZc
        validtriangle = false
    elseif iZa == iZb == iZc
        validtriangle = false
    elseif rZa == rZb
        iZa == iZc || iZb == iZc || (validtriangle = false)
    elseif rZa == rZc
        iZa == iZb || iZb == iZc || (validtriangle = false)
    elseif rZb == rZc
        iZa == iZb || iZa == iZc || (validtriangle = false)
    else
        validtriangle = false
    end
    validtriangle || throw(ArgumentError("`Za`, `Zb`, `Zc` do not define a grid-aligned right triangle"))

    iZa == iZb || ((Zb, Zc) = (Zc, Zb))
    rZb > rZa || ((Za, Zb) = (Zb, Za))
    ==#

    # Determine `dx` and `dy`
    X = rZb - rZa
    Y = abs(iZa - iZc)

    n = ceil(Int, Y/Δr)
    dy = Y/n

    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m

    # precalculate
    slope = Y/X

    v = Vector{complex(typeof(dx))}()
    for j = 0:m  # col
        for i = 0:n  # row (we're traversing down column)

            x = rZa + dx*j
            y = iZc + dy*i

            if y >= (iZc + slope*(x - rZa))
                push!(v, complex(x, y))
            end
        end
    end

    return v
end

"""
    eiwgdomain(Zb, Ze, d, Δr)

Generate a vector of complex coordinates as a rectangulardomain with Δr on the left
and a uppertriangularrectdomain with Δr/5 on the right. The transition occurs at `d` in real.
"""
function eiwgdomain(Zb, Ze, d, Δr)
    if real(Zb) > d
        tricoords = uppertriangularrectdomain(Zb, Ze, Δr/5)
    else
        tricoords = uppertriangularrectdomain(complex(d, imag(Zb)), Ze, Δr/5)
        rectcoords = rectangulardomain(Zb, complex(d, imag(Ze)), Δr)
    end

    origcoords = vcat(tricoords, rectcoords)
    unique!(origcoords)

    return origcoords
end

function uppertriangularrectdomain(Zb::Complex, Ze::Complex, dx, dy)
    rZb, iZb = reim(Zb)
    rZe, iZe = reim(Ze)

    # Determine `dx` and `dy`
    X = rZe - rZb
    Y = abs(iZe - iZb)

    n = ceil(Int, Y/dy)
    dy = Y/n

    m = ceil(Int, X/sqrt(dx^2 - dy^2/4))
    dx = X/m

    slope = 1  # 45° angle

    v = Vector{complex(typeof(X))}()
    for j = 0:m  # col
        x = rZb + dx*j

        for i = 0:n  # row (we're traversing down column)
            y = iZb + dy*i

            if x <= (rZe - (iZe - y)/slope)
                push!(v, complex(x, y))
            end
        end
    end

    return v
end

function uppertriangularrectdomain(Zb::Complex, Ze::Complex, Δr)
    rZb, iZb = reim(Zb)
    rZe, iZe = reim(Ze)

    # Determine `dx` and `dy`
    X = rZe - rZb
    Y = abs(iZe - iZb)

    n = ceil(Int, Y/Δr)
    dy = Y/n

    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m

    slope = 1  # 45° angle

    v = Vector{complex(typeof(X))}()
    for j = 0:m  # col
        x = rZb + dx*j

        for i = 0:n  # row (we're traversing down column)
            y = iZb + dy*i

            if x <= (rZe - (iZe - y)/slope)
                push!(v, complex(x, y))
            end
        end
    end

    return v
end

# TODO: Finish this function
function targetdomain(modes::Vector{T}, delta::T, Δr) where {T<:Complex}
    v = Vector{T}(undef)
    for i in eachindex(modes)
        zb = v[i] - delta
        ze = v[i] + delta
        append!(rectangulardomain(zb, ze, Δr))
    end
    return v
end
