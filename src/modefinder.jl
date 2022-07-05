#==
Functions related to identifying resonant modes (eigenangles) within the earth-ionosphere
waveguide.

These functions use Budden's model, which represents the ground and ionosphere as sharply
reflecting boundaries. The middle of the waveguide is filled with a fictitious medium with a
refractive index that mimicks the propagation of the radio wave through free space over
curved earth.
==#

@doc raw"""
    PhysicalModeEquation{W<:HomogeneousWaveguide} <: ModeEquation

Parameters for solving the physical mode equation ``\det(Rg*R - I)``.

Fields:

    - θ::ComplexF64
    - frequency::Float64
    - waveguide::W
"""
struct PhysicalModeEquation{W<:HomogeneousWaveguide} <: ModeEquation
    θ::ComplexF64
    frequency::Float64
    waveguide::W
end

"""
    PhysicalModeEquation(f, w::HomogeneousWaveguide)

Create a `PhysicalModeEquation` struct with `θ = 0.0+0.0im`.

See also: [`setea`](@ref)
"""
PhysicalModeEquation(f, w::HomogeneousWaveguide) = PhysicalModeEquation(0.0+0.0im, f, w)

"""
    setea(θ, modeequation)

Return `modeequation` with angle `θ`.
"""
setea(θ, modeequation::PhysicalModeEquation) =
    PhysicalModeEquation(θ, modeequation.frequency, modeequation.waveguide)


##########
# Reflection coefficients
##########

"""
    wmatrix(θ, T)

Compute the four submatrix elements of `W` used in the equation ``dR/dz`` from the
ionosphere with `T` matrix returned as a tuple `(W₁₁, W₂₁, W₁₂, W₂₂)`.

Following Budden's formalism for the reflection matrix of a plane wave obliquely incident on
the ionosphere [Budden1955a], the wave below the ionosphere can be resolved into
upgoing and downgoing waves of elliptical polarization, each of whose components are
themselves resolved into a component with the electric field in the plane of propagation and
a component perpendicular to the plane of propagation. The total field can be written in
matrix form as ``e = Lf`` where ``L`` is a 4×4 matrix that simply selects and specifies the
incident angle of the components and ``f`` is a column matrix of the complex amplitudes of
the component waves. By inversion, ``f = L⁻¹e`` and its derivative with respect to height
``z`` is ``f′ = -iL⁻¹TLf = -½iWf``. Then ``W = 2L⁻¹TL`` describes the change in amplitude of
the upgoing and downgoing component waves.

``W`` is also known as ``S`` in many texts.

# References

[Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing
    reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227,
    no. 1171, pp. 516–537, Feb. 1955.
"""
function wmatrix(θ, T)
    C = cos(θ)
    Cinv = 1/C  # == sec(θ)

    # Precompute
    T12Cinv = T[1,2]*Cinv
    T14Cinv = T[1,4]*Cinv
    T32Cinv = T[3,2]*Cinv
    T34Cinv = T[3,4]*Cinv
    CT41 = C*T[4,1]

    #==
    T = SMatrix{4,4}(T11, 0, T31, T41,
                     T12, 0, T32, T42,
                     0,   1,   0,   0,
                     T14, 0, T34, T44)
    L = SMatrix{4,4}(C,  0,  0, 1,
                     0, -1, -C, 0,
                    -C,  0,  0, 1,
                     0, -1,  C, 0)
    # 2*inv(L) = Linv2
    L2inv = SMatrix{4,4}(Cinv,  0, -Cinv,    0,
                         0,    -1,     0,   -1,
                         0, -Cinv,     0, Cinv,
                         1,     0,     1,    0)
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
    dwmatrix(θ, T, dT)

Compute the four submatrix elements of ``dW/dθ`` returned as the tuple
`(dW₁₁, dW₂₁, dW₁₂, dW₂₂)` from the ionosphere with `T` matrix and its derivative with
respect to ``θ``, `dT`.
"""
function dwmatrix(θ, T, dT)
    S, C = sincos(θ)
    Cinv = 1/C
    C²inv = Cinv^2

    dC = -S
    dCinv = S*C²inv

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

Compute the differential of the reflection matrix `R`, ``dR/dz``, at height `z`.

Following the Budden formalism for the reflection of an (obliquely) incident plane wave from
a horizontally stratified ionosphere [Budden1955a], the differential of the
reflection matrix `R` with height `z` can be described by
```math
dR/dz = k/(2i)⋅(W₂₁ + W₂₂R - RW₁₁ - RW₁₂R)
```
Integrating ``dR/dz`` downwards through the ionosphere gives the reflection matrix ``R`` for
the ionosphere as if it were a sharp boundary at the stopping level with free space below.

# References

[Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing
    reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227,
    no. 1171, pp. 516–537, Feb. 1955.
"""
function dRdz(R, modeequation, z; params=LMPParams())
    @unpack θ, frequency = modeequation

    k = wavenumber(frequency)

    M = susceptibility(z, modeequation; params)
    T = tmatrix(θ, M)
    W11, W21, W12, W22 = wmatrix(θ, T)

    # the factor k/(2i) isn't explicitly in [Budden1955a] because of his change of variable
    # ``s = kz``
    return k/2im*(W21 + W22*R - R*W11 - R*W12*R)
end

"""
    dRdθdz(RdRdθ, p, z)

Compute the differential ``dR/dθ/dz`` at height `z` returned as an `SMatrix{4,2}` with
``dR/dz`` in the first 2 rows and ``dR/dθ/dz`` in the bottom 2 rows.

`p` is a tuple containing instances `(PhysicalModeEquation(), LMPParams())`.
"""
function dRdθdz(RdRdθ, p, z)
    modeequation, params = p
    @unpack θ, frequency = modeequation

    k = wavenumber(frequency)

    M = susceptibility(z, modeequation; params)
    T = tmatrix(θ, M)
    dT = dtmatrix(θ, M)
    W11, W21, W12, W22 = wmatrix(θ, T)
    dW11, dW21, dW12, dW22 = dwmatrix(θ, T, dT)

    R = RdRdθ[SVector(1,2),:]
    dRdθ = RdRdθ[SVector(3,4),:]

    dz = k/2im*(W21 + W22*R - R*W11 - R*W12*R)
    dθdz = k/2im*(dW21 + dW22*R + W22*dRdθ - (dRdθ*W11 + R*dW11) -
                    (dRdθ*W12*R + R*dW12*R + R*W12*dRdθ))

    return vcat(dz, dθdz)
end

"""
    integratedreflection(modeequation::PhysicalModeEquation; params=LMPParams())

Integrate ``dR/dz`` downward through the ionosphere described by `modeequation` from
`params.topheight`, returning the ionosphere reflection coefficient `R` at the ground.

`params.integrationparams` are passed to `DifferentialEquations.jl`.
"""
function integratedreflection(modeequation::PhysicalModeEquation; params=LMPParams())

    @unpack topheight, integrationparams = params
    @unpack tolerance, solver, dt, force_dtmin, maxiters = integrationparams

    Mtop = susceptibility(topheight, modeequation; params)
    Rtop = bookerreflection(modeequation.θ, Mtop)

    prob = ODEProblem{false}((R,p,z)->dRdz(R,p,z; params), Rtop, (topheight, BOTTOMHEIGHT), modeequation)

    # WARNING: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, solver; abstol=tolerance, reltol=tolerance,
                force_dtmin=force_dtmin, dt=dt, maxiters=maxiters,
                save_on=false, save_start=false, save_end=true)

    R = sol[end]

    return R
end

"""
    integratedreflection(modeequation::PhysicalModeEquation, ::Dθ; params=LMPParams())

Compute ``R`` and ``dR/dθ`` as an `SMatrix{4,2}` with ``R`` in rows (1, 2) and ``dR/dθ`` in
rows (3, 4).

The `params.integrationparams.tolerance` is hardcoded to `1e-10` in this version of the
function.
"""
function integratedreflection(modeequation::PhysicalModeEquation, ::Dθ; params=LMPParams())
    @unpack topheight, integrationparams = params
    @unpack solver, dt, force_dtmin, maxiters = integrationparams

    # Tolerance is overridden for this `::Dθ` form.
    # Using an identical accuracy appears to result in relatively less accurate solutions
    # compared to the non-Dθ form.
    tolerance = 1e-10

    Mtop = susceptibility(topheight, modeequation; params)
    Rtop, dRdθtop = bookerreflection(modeequation.θ, Mtop, Dθ())
    RdRdθtop = vcat(Rtop, dRdθtop)

    prob = ODEProblem{false}(dRdθdz, RdRdθtop, (topheight, BOTTOMHEIGHT),
                             (modeequation, params))

    # WARNING: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, solver; abstol=tolerance, reltol=tolerance,
                force_dtmin=force_dtmin, dt=dt, maxiters=maxiters,
                save_on=false, save_start=false, save_end=true)

    RdR = sol[end]

    return RdR
end

##########
# Ground reflection coefficient matrix
##########

"""
    fresnelreflection(θ, ground::Ground, frequency)
    fresnelreflection(m::PhysicalModeEquation)

Compute the Fresnel reflection coefficient matrix for the ground-freespace interface at the
ground for a wave `frequency` in Hertz.
"""
fresnelreflection

function fresnelreflection(θ, ground::Ground, frequency)
    S, C = sincos(θ)
    S² = S^2
    ω = angular(frequency)

    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*E0))

    CNg² = C*Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    Rg11 = (CNg² - sqrtNg²mS²)/(CNg² + sqrtNg²mS²)
    Rg22 = (C - sqrtNg²mS²)/(C + sqrtNg²mS²)

    Rg = SDiagonal(Rg11, Rg22)

    return Rg
end

fresnelreflection(m::PhysicalModeEquation) =
    fresnelreflection(m.θ, m.waveguide.ground, m.frequency)

"""
    fresnelreflection(θ, ground::Ground, frequency, ::Dθ)
    fresnelreflection(m::PhysicalModeEquation, ::Dθ)

Compute the Fresnel reflection coefficient matrix for the ground as well as its derivative
with respect to ``θ`` returned as the tuple `(Rg, dRg)`.
"""
function fresnelreflection(θ, ground::Ground, frequency, ::Dθ)
    S, C = sincos(θ)
    S² = S^2
    S2 = 2*S
    ω = angular(frequency)

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
    fresnelreflection(m.θ, m.waveguide.ground, m.frequency, Dθ())

##########
# Identify EigenAngles
##########

"""
    modalequation(R, Rg)

Compute the determinental mode equation ``det(Rg R - I)`` given reflection coefficients `R`
and `Rg`.

A propagating waveguide mode requires that a wave, having reflected from the ionosphere and
then the ground, must be identical with the original upgoing wave. This criteria is met at
roots of the mode equation [Budden1962].

# References

[Budden1962]: K. G. Budden and N. F. Mott, “The influence of the earth’s magnetic field on
    radio propagation by wave-guide modes,” Proceedings of the Royal Society of London.
    Series A. Mathematical and Physical Sciences, vol. 265, no. 1323, pp. 538–553,
    Feb. 1962.
"""
function modalequation(R, Rg)
    return det(Rg*R - I)
end

"""
    dmodalequation(R, dR, Rg, dRg)

Compute the derivative of the determinantal mode equation with respect to ``θ``.
"""
function dmodalequation(R, dR, Rg, dRg)
    # See e.g. https://folk.ntnu.no/hanche/notes/diffdet/diffdet.pdf
    A = Rg*R - I
    dA = dRg*R + Rg*dR
    return det(A)*tr(A\dA)
end

"""
    solvemodalequation(modeequation::PhysicalModeEquation; params=LMPParams())

Compute the ionosphere and ground reflection coefficients and return the value of the
determinental modal equation associated with `modeequation`.

See also: [`solvedmodalequation`](@ref)
"""
function solvemodalequation(modeequation::PhysicalModeEquation; params=LMPParams())
    R = integratedreflection(modeequation; params)
    Rg = fresnelreflection(modeequation)

    f = modalequation(R, Rg)
    return f
end

"""
    solvemodalequation(θ, modeequation::PhysicalModeEquation; params=LMPParams())

Set `θ` for `modeequation` and then solve the modal equation.
"""
function solvemodalequation(θ, modeequation::PhysicalModeEquation; params=LMPParams())
    modeequation = setea(θ, modeequation)
    solvemodalequation(modeequation; params)
end

"""
    solvedmodalequation(modeequation::PhysicalModeEquation; params=LMPParams())

Compute the derivative of the modal equation with respect to ``θ`` returned as the tuple
`(dF, R, Rg)` for the ionosphere and ground reflection coefficients.
"""
function solvedmodalequation(modeequation::PhysicalModeEquation; params=LMPParams())
    RdR = integratedreflection(modeequation, Dθ(); params)
    R = RdR[SVector(1,2),:]
    dR = RdR[SVector(3,4),:]

    Rg, dRg = fresnelreflection(modeequation, Dθ())

    dF = dmodalequation(R, dR, Rg, dRg)
    return dF, R, Rg
end

"""
    solvedmodalequation(θ, modeequation::PhysicalModeEquation; params=LMPParams())

Set `θ` for `modeequation` and then solve the derivative of the mode equation with respect
to `θ`.
"""
function solvedmodalequation(θ, modeequation::PhysicalModeEquation; params=LMPParams())
    modeequation = setea(θ, modeequation)
    solvedmodalequation(modeequation; params)
end

"""
    secantmethod(θ₀, modeequation; params=LMPParams(), tolerance=1e-10, maxiter=200)

Apply secant method to refine eigenangle `θ₀`. `tolerance` is the largest acceptable change
in `θ` between iterations and `maxiter` is the maximum number of iterations before returning
the refined eigenangle. 
"""
function secantmethod(θ₀, modeequation; params=LMPParams())
    @unpack refineeigenangles_tolerance, refineeigenangles_maxiter = params
    tolerance = refineeigenangles_tolerance
    maxiter = refineeigenangles_maxiter

    tolerance² = tolerance^2  # `abs2` is faster than `abs` for complex types
    t = params.grpfparams.tolerance
    δθ = complex(t, t)
    θ₁ = θ₀ + δθ
    iter = 0
    let θ₂
        while abs2(θ₁ - θ₀) > tolerance² && iter < maxiter
            f₀ = solvemodalequation(θ₀, modeequation; params)
            f₁ = solvemodalequation(θ₁, modeequation; params)
            θ₂ = θ₁ - f₁*(θ₁ - θ₀)/(f₁ - f₀)
            θ₀, θ₁ = θ₁, θ₂
            iter += 1
        end
        iter >= maxiter && @warn "maxiter reached"
        return θ₂
    end
end

"""
    findmodes(modeequation::ModeEquation, mesh=nothing; params=LMPParams())

Find eigenangles associated with `modeequation.waveguide` within the domain of
`mesh`.

`mesh` should be an array of complex numbers that make up the original grid over which
the GRPF algorithm searches for roots of `modeequation`. If `mesh === nothing`,
it is computed with [`defaultmesh`](@ref).

There is a check for redundant modes that requires modes to be separated by at least
1 orders of magnitude greater than `grpfparams.tolerance` in real and/or imaginary
component. For example, if `grpfparams.tolerance = 1e-4`, then either the real or imaginary
component of each mode must be separated by at least 1e-3 radians from every other mode.
"""
function findmodes(modeequation::ModeEquation, mesh=nothing; params=LMPParams())
    @unpack grpfparams, refineeigenangles = params

    if isnothing(mesh)
        mesh = defaultmesh(modeequation.frequency)
    end

    # WARNING: If tolerance of mode finder is much less than the R integration tolerance, it
    # is possible multiple identical modes will be identified. Checks for valid and
    # redundant modes help ensure valid eigenangles are returned from this function.

    roots, _ = grpf(θ->solvemodalequation(θ, modeequation; params), mesh, grpfparams)

    # Scale tolerance for filtering
    # if tolerance is 1e-4, this rounds to 3 decimal places
    ndigits = round(Int, abs(log10(grpfparams.tolerance)+1), RoundDown)

    # Remove any redundant modes
    sort!(roots; by=reim, rev=true)
    unique!(z->round(z; digits=ndigits), roots)

    if refineeigenangles
        # use higher tolerance for integrationparams
        for i in eachindex(roots)
            θ = secantmethod(roots[i], modeequation;
                params=LMPParams(params; integrationparams=IntegrationParams(tolerance=1e-8)))
            roots[i] = θ
        end
    end

    return roots
end

#==
Mesh grids for `GRPF`
==#

"""
    defaultmesh(frequency; rmin=deg2rad(30.0), imin=deg2rad(-10.0),
        Δr_coarse=deg2rad(0.5), Δr_fine=deg2rad(0.1),
        rtransition=deg2rad(75.0), itransition=deg2rad(-1.5))

Generate vector of complex coordinates to be used by GRPF in the search for
waveguide modes.

`rmin` is the lower bound of the real axis and `imin` is the lower bound of the imaginary
axis.

At frequencies above 12 kHz the mesh spacing in the upper right corner of the domain
with real values above `rtransition` and imaginary values above `itransition` is
`Δr_fine` and is `Δr_coarse` everywhere else.

At frequencies below 12 kHz the mesh spacing is always `Δr_coarse`.

The lower right diagonal of the lower right quadrant of the complex plane is excluded
from the mesh.

See also: [`findmodes`](@ref)
"""
function defaultmesh(frequency;
    rmin=deg2rad(30.0), imin=deg2rad(-10.0),
    Δr_coarse=deg2rad(0.5), Δr_fine=deg2rad(0.1),
    rtransition=deg2rad(75.0), itransition=deg2rad(-1.5))

    # TODO: get a better idea of frequency transition
    if frequency > 12000
        zbl_coarse = complex(rmin, imin)
        ztr_coarse = complex(deg2rad(89.9), 0.0)

        mesh = trianglemesh(zbl_coarse, ztr_coarse, Δr_coarse)

        filter!(z->(real(z) < rtransition || imag(z) < itransition), mesh)

        zbl_fine = complex(rtransition, itransition)
        ztr_fine = complex(deg2rad(89.9), 0.0)

        append!(mesh, trianglemesh(zbl_fine, ztr_fine, Δr_fine))
    else
        zbl = complex(rmin, imin)
        ztr = complex(deg2rad(89.9), 0.0)

        mesh = trianglemesh(zbl, ztr, Δr_coarse)
    end

    return mesh
end

"""
    trianglemesh(zbl, ztr, Δr)

Generate initial mesh node coordinates for a right triangle domain from complex
coordinate `zbl` in the bottom left and `ztr` in the top right with initial mesh step `Δr`.

`zbl` and `ztr` are located on the complex plane at the locations marked in the
diagram below:

```
       im
       |
-re ---|------ ztr
       |     /
       |    /
       |   /
       |  /
       zbl
```
"""
function trianglemesh(zbl, ztr, Δr)
    rzbl, izbl = reim(zbl)
    rztr, iztr = reim(ztr)

    X = rztr - rzbl
    Y = iztr - izbl

    n = ceil(Int, Y/Δr)
    dy = Y/n

    ## dx = sqrt(Δr² - (dy/2)²), solved for equilateral triangle
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m
    half_dx = dx/2

    slope = 1  # 45° angle

    T = promote_type(ComplexF64, typeof(zbl), typeof(ztr), typeof(Δr))
    mesh = Vector{T}()

    shift = false  # we will displace every other line by dx/2
    for j = 0:n
        y = izbl + dy*j

        for i = 0:m
            x = rzbl + dx*i

            if shift && i == 0
                continue  # otherwise, we shift out of left bound
            elseif shift
                x -= half_dx
            end

            if i == m
                shift = !shift
            end

            ## NEW: check if `x, y` is in upper left triangle
            if y >= slope*x - π/2
                push!(mesh, complex(x, y))
            end
        end
    end

    return mesh
end
