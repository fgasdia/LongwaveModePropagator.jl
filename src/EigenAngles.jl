# EigenAngle.jl
# 
# Plane-wave propagation direction in angle `θ` from the vertical in radians (because df/dθ
# are in radians).

"""
    isdetached(θ, frequency; params=LMPParams())

Return `true` if angle `θ` is likely an earth detached (whispering gallery) mode according to the
criteria in [Pappert1981] eq. 1 with the additional criteria that the frequency be
above 100 kHz.

# References

[Pappert1981]: R. A. Pappert, “LF daytime earth ionosphere waveguide calculations,” Naval
    Ocean Systems Center, San Diego, CA, NOSC/TR-647, Jan. 1981.
    [Online]. Available: https://apps.dtic.mil/docs/citations/ADA096098.
"""
function isdetached(θ, frequency; params=LMPParams())
    @unpack earthradius, curvatureheight = params

    C = cos(θ)
    C² = C^2
    k = wavenumber(frequency)
    α = 2/earthradius

    return frequency > 100e3 && real(2im/3*(k/α)*(C²- α*curvatureheight)^(3/2)) > 12.4
end

"""
    referencetoground(ea; params=LMPParams())

Reference eigenangle `ea` from `params.curvatureheight` to ground (``z = 0``).
"""
function referencetoground(ea; params=LMPParams())
    # see, e.g. PS71 pg 11
    @unpack earthradius, curvatureheight = params
    return asin(sin(ea)/sqrt(1 - 2/earthradius*curvatureheight))
end

"""
    attenuation(ea, frequency)

Compute attenuation of eigenangle `ea` at the ground for a wave `frequency` in Hertz.

This function internally references `ea` to the ground.
"""
function attenuation(ea, frequency; params=LMPParams())
    ea₀ = referencetoground(ea; params)
    S₀ = sin(ea₀)
    neper2dB = 20log10(exp(1))  # 1 Np ≈ 8.685 dB
    return -neper2dB*wavenumber(frequency)*imag(S₀)*1e6
end

"""
    phasevelocity(ea; params=LMPParams())

Compute the relative phase velocity ``v/c`` associated with the eigenangle `ea`.

This function internally references `θ` to the ground.
"""
function phasevelocity(ea; params=LMPParams())
    ea₀ = referencetoground(ea; params)
    S₀ = sin(ea₀)
    return 1/real(S₀)
end
