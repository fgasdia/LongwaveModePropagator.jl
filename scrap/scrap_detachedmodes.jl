
"""
    detachedreflection(ea, frequency::Frequency; params=LMPParams())

Compute the ground reflection coefficient for earth detached (whispering gallery) modes.

# References

[^Pappert1981]: R. A. Pappert, “LF daytime earth ionosphere waveguide calculations,” Naval
    Ocean Systems Center, San Diego, CA, NOSC/TR-647, Jan. 1981.
    [Online]. Available: https://apps.dtic.mil/docs/citations/ADA096098.
"""
function detachedreflection(ea, frequency::Frequency; params=LMPParams())
    ea = EigenAngle(ea)

    @unpack earthradius, curvatureheight = params
    C, C² = ea.cosθ, ea.cos²θ
    k = frequency.k

    α = 2/earthradius

    q₀ = pow23(k/α)*(C²- α*curvatureheight)
    h₁0, h₂0, dh₁0, dh₂0 = modifiedhankel(q₀)
    h = h₂0 - exp(4im*π/3)*h₁0
    dh = dh₂0 - exp(4im*π/3)*dh₁0

    n₀² = 1 - α*curvatureheight

    tmp = cbrt(α/k)*dh
    F = 1im/n₀²*((α/k)*h/2 + tmp)  # = 1im/n₀²*((α/k)*h/2 + cbrt(α/k)*dh)

    Rg11 = (C*h + F)/(C*h - F)
    Rg22 = (C*h + 1im*tmp)/(C*h - 1im*tmp)

    Rg = SDiagonal(Rg11, Rg22)

    return Rg
end

"`detachedreflection(m::PhysicalModeEquation; params=LMPParams())`"
detachedreflection(m::PhysicalModeEquation; params=LMPParams()) =
    detachedreflection(m.ea, m.frequency, params=params)

"""
    detachedreflection(ea, frequency::Frequency, ::Dθ; params=LMPParams())

Compute the reflection coefficient matrix for the ground for earth detached modes
[^Pappert1981] as well as its derivative with respect to ``θ`` returned as the tuple
`(Rg, dRg)`.

# References

[^Pappert1981]: R. A. Pappert, “LF daytime earth ionosphere waveguide calculations,” Naval
    Ocean Systems Center, San Diego, CA, NOSC/TR-647, Jan. 1981.
    [Online]. Available: https://apps.dtic.mil/docs/citations/ADA096098.
"""
function detachedreflection(ea, frequency::Frequency, ::Dθ; params=LMPParams())
    ea = EigenAngle(ea)
    Rg = detachedreflection(ea, frequency, params=params)

    # TODO: derivative through modifiedhankel
    dRg = FiniteDiff.finite_difference_derivative(θ->detachedreflection(θ, frequency), ea.θ,
                                                  Val{:central}, ComplexF64, Rg)

    return Rg, dRg
end

"`detachedreflection(m::PhysicalModeEquation, ::Dθ; params=LMPParams())`"
detachedreflection(m::PhysicalModeEquation, ::Dθ; params=LMPParams()) =
    detachedreflection(m.ea, m.frequency, Dθ(), params=params)
