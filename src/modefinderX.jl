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

struct ModifiedModeEquation{W<:HomogeneousWaveguide,F} <: ModeEquation
    frequency::Frequency
    waveguide::W
    Mfcn::F
end

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
function dXdz(X, params, z)
    @unpack ea, frequency, Mfcn = params

    k = frequency.k
    C = ea.cosθ
    C² = ea.cos²θ

    M = Mfcn(z)
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

function integratedreflectionX(ea::EigenAngle, frequency::Frequency,
    waveguide::HomogeneousWaveguide, Mfcn::F) where F

    @unpack bfield, species = waveguide
    Mtop = susceptibility(TOPHEIGHT, frequency, bfield, species)
    Rtop = sharpboundaryreflection(ea, Mtop)
    Xtop = R2X(ea, Rtop)

    dzparams = DZParams(ea, frequency, Mfcn)
    prob = ODEProblem{false}(dXdz, Xtop, (TOPHEIGHT, BOTTOMHEIGHT), dzparams)

    #==
    | tolerance | method |
    |-----------|--------|
    | 1e-8      | Vern8  | slightly higher memory, much faster than Vern6, Vern7 twice as slow as 1e-6
    | 1e-6      | Vern8 |
    ==#

    # NOTE: When save_on=false, don't try interpolating the solution!
    sol = solve(prob, Vern8(), abstol=1e-8, reltol=1e-8,
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

function solvemodalequation(ea::EigenAngle, modeequation::ModifiedModeEquation)
    @unpack frequency, waveguide, Mfcn = modeequation
    X = integratedreflectionX(ea, frequency, waveguide, Mfcn)
    n, d = fresnelreflectionX(ea, waveguide.ground, frequency)

    f = modalequationX(X, n, d)
    return f
end
