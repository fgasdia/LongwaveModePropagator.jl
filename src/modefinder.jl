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

    - ea::EigenAngle
    - frequency::Frequency
    - waveguide::W
"""
struct PhysicalModeEquation{W<:HomogeneousWaveguide} <: ModeEquation
    ea::EigenAngle
    frequency::Frequency
    waveguide::W
end

"""
    PhysicalModeEquation(f::Frequency, w::HomogeneousWaveguide)

Create a `PhysicalModeEquation` struct with `ea = complex(0.0)`.

See also: [`setea`](@ref)
"""
PhysicalModeEquation(f::Frequency, w::HomogeneousWaveguide) =
    PhysicalModeEquation(EigenAngle(complex(0.0)), f, w)

"""
    setea(ea, modeequation)

Return `modeequation` with eigenangle `ea`.

`ea` will be converted to an `EigenAngle` if necessary.
"""
setea(ea, modeequation::PhysicalModeEquation) =
    PhysicalModeEquation(EigenAngle(ea), modeequation.frequency, modeequation.waveguide)

"""
    isroot(x; atol=1e-2)

Return `true` if `x` is approximately equal to 0 with the absolute tolerance `atol`.
"""
isroot(x; atol=1e-2, kws...) = isapprox(x, 0, atol=atol, kws...)

"""
    filterroots!(roots, frequency, waveguide; atol=0.1)
    filterroots!(roots, modeequation; atol=0.1)

Remove elements from `roots` if they are not valid roots of the physical modal equation.

Although the default `atol` for this function is larger than for `isroot`(@ref), the
modal equation is very sensitive to θ, and thus is much larger than 0.1 if θ is off
from an eigenangle even slightly.
"""
filterroots!

function filterroots!(roots, frequency, waveguide; atol=0.1)
    modeequation = PhysicalModeEquation(frequency, waveguide)
    return filterroots!(roots, modeequation, atol=atol)
end

function filterroots!(roots, modeequation::PhysicalModeEquation; atol=0.1)
    i = 1
    while i <= length(roots)
        modeequation = setea(EigenAngle(roots[i]), modeequation)
        f = solvemodalequation(modeequation)
        isroot(f, atol=atol) ? (i += 1) : deleteat!(roots, i)
    end
    return roots
end

##########
# Reflection coefficients
##########

"""
    wmatrix(ea::EigenAngle, T)

Compute the four submatrix elements of `W` used in the equation ``dR/dz`` from the
ionosphere with `T` matrix returned as a tuple `(W₁₁, W₂₁, W₁₂, W₂₂)`.

Following Budden's formalism for the reflection matrix of a plane wave obliquely incident on
the ionosphere [[Budden1955a](@cite)], the wave below the ionosphere can be resolved into
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
function wmatrix(ea::EigenAngle, T)
    C, Cinv = ea.cosθ, ea.secθ

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
    dwmatrix(ea::EigenAngle, T, dT)

Compute the four submatrix elements of ``dW/dθ`` returned as the tuple
`(dW₁₁, dW₂₁, dW₁₂, dW₂₂)` from the ionosphere with `T` matrix and its derivative with
respect to ``θ``, `dT`.
"""
function dwmatrix(ea::EigenAngle, T, dT)
    C, S, C², Cinv = ea.cosθ, ea.sinθ, ea.cos²θ, ea.secθ
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
    dRdz(R, p, z)

Compute the differential of the reflection matrix `R`, ``dR/dz``, at height `z`. `p` is a
tuple containing instances `(PhysicalModeEquation(), LMPParams())`.

Following the Budden formalism for the reflection of an (obliquely) incident plane wave from
a horizontally stratified ionosphere [[Budden1955a](@cite)], the differential of the
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
function dRdz(R, p, z)
    modeequation, params = p
    @unpack ea, frequency = modeequation

    k = frequency.k

    M = susceptibility(z, modeequation, params=params)
    T = tmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)

    # the factor k/(2i) isn't explicitly in [^Budden1955a] because of his change of variable
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
    @unpack ea, frequency = modeequation

    k = frequency.k

    M = susceptibility(z, modeequation, params=params)
    T = tmatrix(ea, M)
    dT = dtmatrix(ea, M)
    W11, W21, W12, W22 = wmatrix(ea, T)
    dW11, dW21, dW12, dW22 = dwmatrix(ea, T, dT)

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

    Mtop = susceptibility(topheight, modeequation, params=params)
    Rtop = bookerreflection(modeequation.ea, Mtop)

    prob = ODEProblem{false}(dRdz, Rtop, (topheight, BOTTOMHEIGHT), (modeequation, params))

    # WARNING: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, solver, abstol=tolerance, reltol=tolerance,
                force_dtmin=force_dtmin, dt=dt, maxiters=maxiters,
                save_on=false, save_start=false, save_end=true)

    R = sol[end]

    return R
end

"""
    integratedreflection(modeequation::PhysicalModeEquation, ::Dθ; params=LMPParams())

Compute ``R`` and ``dR/dθ`` as an `SMatrix{4,2}` with ``R`` in rows (1, 2) and ``dR/dθ`` in
rows (3, 4).
"""
function integratedreflection(modeequation::PhysicalModeEquation, ::Dθ; params=LMPParams())
    @unpack topheight, integrationparams = params
    @unpack tolerance, solver, dt, force_dtmin, maxiters = integrationparams

    Mtop = susceptibility(topheight, modeequation, params=params)
    Rtop, dRdθtop = bookerreflection(modeequation.ea, Mtop, Dθ())
    RdRdθtop = vcat(Rtop, dRdθtop)

    prob = ODEProblem{false}(dRdθdz, RdRdθtop, (topheight, BOTTOMHEIGHT),
                             (modeequation, params))

    # WARNING: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, solver, abstol=tolerance, reltol=tolerance,
                force_dtmin=force_dtmin, dt=dt, maxiters=maxiters,
                save_on=false, save_start=false, save_end=true)

    RdR = sol[end]

    return RdR
end

##########
# Ground reflection coefficient matrix
##########

"""
    fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency)
    fresnelreflection(m::PhysicalModeEquation)

Compute the Fresnel reflection coefficient matrix for the ground-freespace interface at the
ground.
"""
fresnelreflection

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

"""
    fresnelreflection(ea::EigenAngle, ground::Ground, frequency::Frequency, ::Dθ)
    fresnelreflection(m::PhysicalModeEquation, ::Dθ)

Compute the Fresnel reflection coefficient matrix for the ground as well as its derivative
with respect to ``θ`` returned as the tuple `(Rg, dRg)`.
"""
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

Compute the determinental mode equation ``det(Rg R - I)`` given reflection coefficients `R`
and `Rg`.

A propagating waveguide mode requires that a wave, having reflected from the ionosphere and
then the ground, must be identical with the original upgoing wave. This criteria is met at
roots of the mode equation [[Budden1962](@cite)].

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
    R = integratedreflection(modeequation, params=params)
    Rg = fresnelreflection(modeequation)

    f = modalequation(R, Rg)
    return f
end

"""
    solvemodalequation(θ, modeequation::PhysicalModeEquation; params=LMPParams())

Set `θ` for `modeequation` and then solve the modal equation.
"""
function solvemodalequation(θ, modeequation::PhysicalModeEquation; params=LMPParams())
    # Convenience function for `grpf`
    modeequation = setea(EigenAngle(θ), modeequation)
    solvemodalequation(modeequation, params=params)
end

"""
    solvedmodalequation(modeequation::PhysicalModeEquation; params=LMPParams())

Compute the derivative of the modal equation with respect to ``θ`` returned as the tuple
`(dF, R, Rg)` for the ionosphere and ground reflection coefficients.
"""
function solvedmodalequation(modeequation::PhysicalModeEquation; params=LMPParams())
    RdR = integratedreflection(modeequation, Dθ(), params=params)
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
    modeequation = setea(EigenAngle(θ), modeequation)
    solvedmodalequation(modeequation, params=params)
end

"""
    findmodes(modeequation::ModeEquation, origcoords=nothing; params=LMPParams())

Find `EigenAngle`s associated with `modeequation.waveguide` within the domain of
`origcoords`.

`origcoords` should be an array of complex numbers that make up the original grid over which
the GRPF algorithm searches for roots of `modeequation`. If `origcoords == nothing`,
it is computed with [`defaultmesh`](@ref).

Roots found by the GRPF algorithm are confirmed to a tolerance of approximately 3 orders
of magnitude greater than `grpfparams.tolerance` in both the real and imaginary component.
For example, if `grpfparams.tolerance = 1e-5`, then the value of the modal equation for each
root must be less than `1e-2`. Typically the values of each component are close to
`grpfparams.tolerance`.

There is also a check for redundant modes that requires modes to be separated by at least
2 orders of magnitude greater than `grpfparams.tolerance` in real and/or imaginary
component. For example, if `grpfparams.tolerance = 1e-5`, then either the real or imaginary
component of each mode must be separated by at least 1e-3 from every other mode.
"""
function findmodes(modeequation::ModeEquation, origcoords=nothing; params=LMPParams())
    @unpack grpfparams = params

    if isnothing(origcoords)
        origcoords = defaultmesh(modeequation.frequency)
    end

    # WARNING: If tolerance of mode finder is much less than the R integration tolerance, it
    # is possible multiple identical modes will be identified. Checks for valid and
    # redundant modes help ensure valid eigenangles are returned from this function.

    roots, poles = grpf(θ->solvemodalequation(θ, modeequation, params=params),
                        origcoords, grpfparams)

    # Scale tolerance for filtering
    # if tolerance is 1e-8, this rounds to 6 decimal places
    ndigits = round(Int, abs(log10(grpfparams.tolerance)+2), RoundDown)

    # Ensure roots are valid solutions to the modal equation
    filtertolerance = exp10(-ndigits+1)  # +3 from the grpfparams.tolerannce
    filterroots!(roots, modeequation, atol=filtertolerance)

    # Remove any redundant modes
    sort!(roots, by=reim, rev=true)
    unique!(z->round(z, digits=ndigits), roots)

    return EigenAngle.(roots)
end

#==
Mesh grids for `GRPF`
==#

"""
    defaultmesh(frequency)

Generate vector of complex coordinates to be used by GRPF in the search for
waveguide modes.

At frequencies above 12 kHz the mesh spacing in the upper right corner of the domain
with real values above 80° and imaginary values above -1° is ``0.15 π/180`` and is
``0.5 π/180`` everywhere else.

The lower right diagonal of the lower right quadrant of the complex plane is excluded
from the mesh.

See also: [`findmodes`](@ref)
"""
function defaultmesh(frequency)

    # TODO: get a better idea of frequency transition
    if frequency > 12000
        zbl_coarse = complex(deg2rad(30.0), deg2rad(-10.0))
        ztr_coarse = complex(deg2rad(89.9), 0.0)
        Δr_coarse = deg2rad(0.5)

        mesh = trianglemesh(zbl_coarse, ztr_coarse, Δr_coarse)

        rtransition = deg2rad(80.0)
        itransition = deg2rad(-1.0)
        filter!(z->(real(z) < rtransition || imag(z) < itransition), mesh)

        zbl_fine = complex(rtransition, itransition)
        ztr_fine = complex(deg2rad(89.9), 0.0)
        Δr_fine = deg2rad(0.15)

        append!(mesh, trianglemesh(zbl_fine, ztr_fine, Δr_fine))
    else
        zbl = complex(deg2rad(30.0), deg2rad(-10.0))
        ztr = complex(deg2rad(89.9), 0.0)
        Δr = deg2rad(0.5)

        mesh = trianglemesh(zbl, ztr, Δr)
    end

    return mesh
end
defaultmesh(f::Frequency) = defaultmesh(f.f)

"""
    trianglemesh(zbl, ztr, Δr)

Generate initial mesh node coordinates for a right triangle domain from complex
coordinate `zbl` in the bottom left and `ztr` in the top right with initial mesh step `Δr`.

`zbl` and `ztr` are located on the complex plane at the locations marked in the
diagram below:

       im
       |
-re ---|------ ztr
       |     /
       |    /
       |   /
       |  /
       zbl

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
