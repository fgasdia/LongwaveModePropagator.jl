#==
Functions related to calculating physical parameters of the ionospheric plasma
(susceptibility) relevant to the propagation of electromagnetic waves.
==#

@doc raw"""
    susceptibility

Compute the ionosphere susceptibility tensor `M` as a `SMatrix{3,3}` using
`species.numberdensity` and `species.collisionfrequency` at `altitude`.

If `params.earthcurvature == true`, `M` includes a first order correction for earth
curvature by means of a fictitious refractive index [^Pappert1967].

The susceptibility matrix is calculated from the constitutive relations presented in
[^Ratcliffe1959]. This includes the effect of earth's magnetic field vector and collisional
damping on electron motion.

The tensor is:
```math
M = -\frac{X}{U(U²-Y²)}
        \begin{pmatrix}
            U² - x²Y²    & -izUY - xyY² & iyUY - xzY² \\
            izUY - xyY²  & U² - y²Y²    & -ixUY - yzY² \\
            -iyUY - xzY² & ixUY - yzY²  & U² - z²Y²
        \end{pmatrix}
```
where ``X = ωₚ²/ω²``, ``Y = |ωₕ/ω|``, ``Z = ν/ω``, and ``U = 1 - iZ``. The earth curvature
correction subtracts ``2/Rₑ*(H - altitude)`` from the diagonal of ``M`` where ``H`` is
`params.curvatureheight`.

# References

[^Pappert1967]: R. A. Pappert, E. E. Gossard, and I. J. Rothmuller, “A numerical
    investigation of classical approximations used in VLF propagation,” Radio Science,
    vol. 2, no. 4, pp. 387–400, Apr. 1967, doi: 10.1002/rds196724387.

[^Ratcliffe1959]: J. A. Ratcliffe, "The magneto-ionic theory & its applications to the
    ionosphere," Cambridge University Press, 1959.
"""
function susceptibility end

"`susceptibility(altitude, frequency, bfield, species; params=LMPParams())`"
function susceptibility(altitude, frequency, bfield, species; params=LMPParams())
    @unpack earthradius, earthcurvature, curvatureheight = params

    B, x, y, z = bfield.B, bfield.dcl, bfield.dcm, bfield.dcn
    x², y², z² = x^2, y^2, z^2
    ω = frequency.ω

    e, m = species.charge, species.mass
    N, nu = species.numberdensity, species.collisionfrequency

    invω = inv(ω)
    invmω = invω/m  # == inv(m*ω)

    # Constitutive relations
    # (see Budden1955a, pg. 517, Budden1988 pg. 39, or [^Ratcliffe1959])
    X = N(altitude)*e^2*invω*invmω/E0
    Y = e*B*invmω  # [^Ratcliffe1959] pg. 182 specifies that sign of `e` is included here
    Z = nu(altitude)*invω
    U = 1 - 1im*Z

    # TODO: Support multiple species, e.g.
    # U²D = zero()
    # Y²D = zero()
    # UYD = zero()
    # for i = 1:length(species)
    #     U²D += U^2*D
    #     UYD += Y*U*D
    #     Y²D += Y^2*D
    # end

    # Precompute variables
    # These are formulated such that application to multiple species is trivial. Direction
    # cosines are applied last
    U² = U^2
    Y² = Y^2

    D = -X/(U*(U² - Y²))

    U²D = U²*D
    Y²D = Y²*D
    UYD = U*Y*D

    # Leverage partial symmetry of M to reduce computations
    izUYD = 1im*z*UYD
    xyY²D = x*y*Y²D
    iyUYD = 1im*y*UYD
    xzY²D = x*z*Y²D
    ixUYD = 1im*x*UYD
    yzY²D = y*z*Y²D

    # Elements of `M`
    M11 = U²D - x²*Y²D
    M21 = izUYD - xyY²D
    M31 = -iyUYD - xzY²D
    M12 = -izUYD - xyY²D
    M22 = U²D - y²*Y²D
    M32 = ixUYD - yzY²D
    M13 = iyUYD - xzY²D
    M23 = -ixUYD - yzY²D
    M33 = U²D - z²*Y²D

    if earthcurvature
        curvaturecorrection = 2/earthradius*(curvatureheight - altitude)

        M11 -= curvaturecorrection
        M22 -= curvaturecorrection
        M33 -= curvaturecorrection
    end

    M = SMatrix{3,3}(M11, M21, M31, M12, M22, M32, M13, M23, M33)

    return M
end

"`susceptibility(altitude, me::ModeEquation; params=LMPParams())`"
susceptibility(altitude, me::ModeEquation; params=LMPParams()) =
    susceptibility(altitude, me.frequency, me.waveguide, params=params)

"`susceptibility(altitude, frequency, w::HomogeneousWaveguide; params=LMPParams())`"
susceptibility(altitude, frequency, w::HomogeneousWaveguide; params=LMPParams()) =
    susceptibility(altitude, frequency, w.bfield, w.species, params=params)
