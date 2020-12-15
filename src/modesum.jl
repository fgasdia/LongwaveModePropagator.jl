#==
Excitation factor, height gain functions, and electric field mode sum
==#

"""
    ExcitationFactor{T,T2}

Constants used in calculating excitation factors and height gains.

# Fields

- `F₁::T`: height gain constant. See [^Pappert1976].
- `F₂::T`
- `F₃::T`
- `F₄::T`
- `h₁0::T`: first modified Hankel function of order 1/3 at the ground.
- `h₂0::T`: second modified Hankel function of order 1/3 at the ground.
- `EyHy::T`: polarization ratio ``Ey/Hy``, derived from reflection coefficients (or ``T``s).
- `Rg::T2`: ground reflection coefficient matrix.

# References

[^Pappert1976]: R. A. Pappert and L. R. Shockey, “Simplified VLF/LF mode conversion program
    with allowance for elevated, arbitrarily oriented electric dipole antennas,” Naval
    Electronics Laboratory Center, San Diego, CA, Interim Report 771, Oct. 1976. [Online].
    Available: http://archive.org/details/DTIC_ADA033412.
"""
struct ExcitationFactor{T,T2}
    F₁::T
    F₂::T
    F₃::T
    F₄::T
    h₁0::T
    h₂0::T
    EyHy::T
    Rg::T2
end

"""
    excitationfactorconstants(ea₀, R, Rg, frequency, ground; params=LMPParams())

Return an `ExcitationFactor` struct used in calculating height-gain functions and excitation
factors where eigenangle `ea₀` is referenced to the ground.

!!! note

    This function assumes that reflection coefficients are referenced to ``d = z = 0``.

# References

[^Pappert1976]: R. A. Pappert and L. R. Shockey, “Simplified VLF/LF mode conversion program
    with allowance for elevated, arbitrarily oriented electric dipole antennas,” Naval
    Electronics Laboratory Center, San Diego, CA, Interim Report 771, Oct. 1976. [Online].
    Available: http://archive.org/details/DTIC_ADA033412.

[^Ferguson1980]: J. A. Ferguson and F. P. Snyder, “Approximate VLF/LF waveguide mode
    conversion model: Computer applications: FASTMC and BUMP,” Naval Ocean Systems Center,
    San Diego, CA, NOSC-TD-400, Nov. 1980. Accessed: May 08, 2017. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA096240.

[^Morfitt1980]: D. G. Morfitt, “‘Simplified’ VLF/LF mode conversion computer programs:
    GRNDMC and ARBNMC,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-514, Jan. 1980.
    Accessed: Jan. 15, 2018. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA082695.
"""
function excitationfactorconstants(ea₀, R, Rg, frequency, ground; params=LMPParams())
    S², C² = ea₀.sin²θ, ea₀.cos²θ
    k, ω = frequency.k, frequency.ω
    ϵᵣ, σ = ground.ϵᵣ, ground.σ

    @unpack earthradius = params

    # Precompute
    α = 2/earthradius
    tmp1 = pow23(α/k)/2  # 1/2*(a/k)^(2/3)

    q₀ = pow23(k/α)*C²  # (a/k)^(-2/3)*C²
    h₁0, h₂0, dh₁0, dh₂0 = modifiedhankel(q₀)

    H₁0 = dh₁0 + tmp1*h₁0
    H₂0 = dh₂0 + tmp1*h₂0

    n₀² = 1  # modified free space index of refraction squared, referenced to ground
    Ng² = complex(ϵᵣ, -σ/(ω*E0))  # ground index of refraction

    # Precompute
    tmp2 = 1im*cbrt(k/α)*sqrt(Ng² - S²)  # i(k/α)^(1/3)*(Ng² - S²)^(1/2)

    F₁ = -H₂0 + (n₀²/Ng²)*tmp2*h₂0
    F₂ = H₁0 - (n₀²/Ng²)*tmp2*h₁0
    F₃ = -dh₂0 + tmp2*h₂0
    F₄ = dh₁0 - tmp2*h₁0

    # ``EyHy = ey/hy``. Also known as `f0fr` or `f`.
    # It is a polarization ratio that adds the proper amount of TE wave when the y component
    # of the magnetic field is normalized to unity at the ground.
    # A principally TE mode will have `1 - R[1,1]*Rg[1,1]` very small and EyHy will be very
    # small, so we use the first equation below. Conversely, a principally TM mode will have
    # `1 - R[2,2]Rg[2,2]` very small and EyHy very large, resulting in the use of the second
    # equation below. [^Ferguson1980] pg. 58 seems to suggest the use of the opposite, but
    # LWPC uses the form used here and this makes sense because there are more working
    # decimal places.
    if abs2(1 - R[1,1]*Rg[1,1]) < abs2(1 - R[2,2]*Rg[2,2])
        # EyHy = T₃/T₁
        EyHy = (1 + Rg[2,2])*R[2,1]*Rg[1,1]/((1 + Rg[1,1])*(1 - R[2,2]*Rg[2,2]))
    else
        # EyHy = T₂/(T₃*T₄)
        EyHy = (1 + Rg[2,2])*(1 - R[1,1]*Rg[1,1])/((1 + Rg[1,1])*R[1,2]*Rg[2,2])
    end

    return ExcitationFactor(F₁, F₂, F₃, F₄, h₁0, h₂0, EyHy, Rg)
end

"""
    excitationfactor(ea, dFdθ, R, Rg, efconstants::ExcitationFactor; params=LMPParams())

Compute excitation factors for the ``Hy`` field at the emitter returned as the tuple
`(λv, λb, λe)` for vertical, broadside, and end-on dipoles.

The excitation factor describes how efficiently the field component can be excited in the
waveguide.

This function most closely follows the approach taken in [^Pappert1983], which makes use of
``T`` rather than ``τ``. From the total ``Hy`` excitation factor (the sum product of the
`λ`s with the antenna orientation terms), the excitation factor for electric fields can
be found as:

- ``λz = -S₀λ``
- ``λx = EyHy⋅λ``
- ``λy = -λ``

!!! note

    This function assumes that `R` and `Rg` are referenced to ``d = z = 0``.

# References

[^Morfitt1980]: D. G. Morfitt, “‘Simplified’ VLF/LF mode conversion computer programs:
    GRNDMC and ARBNMC,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-514, Jan. 1980.
    Accessed: Jan. 15, 2018. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA082695.

[^Pappert1983]: R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low
    Frequency/Very Low Frequency) long path pulse program for antennas of arbitrary
    elevation and orientation,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891,
    Aug. 1983. Accessed: Jul. 04, 2018. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA133876.

[^Pappert1986]: R. A. Pappert and J. A. Ferguson, “VLF/LF mode conversion model calculations
    for air to air transmissions in the earth-ionosphere waveguide,” Radio Sci., vol. 21,
    no. 4, pp. 551–558, Jul. 1986, doi: 10.1029/RS021i004p00551.
"""
function excitationfactor(ea, dFdθ, R, efconstants::ExcitationFactor; params=LMPParams())
    S = ea.sinθ
    sqrtS = sqrt(S)
    S₀ = referencetoground(ea.sinθ, params=params)

    @unpack F₁, F₂, F₃, F₄, h₁0, h₂0, Rg = efconstants

    # Unlike the formulations shown in the references, we scale these excitation factors
    # with `D##` instead of `EyHy` and appropriately don't scale the height gains.
    F₁h₁0 = F₁*h₁0
    F₂h₂0 = F₂*h₂0
    F₃h₁0 = F₃*h₁0
    F₄h₂0 = F₄*h₂0

    D₁₁ = (F₁h₁0 + F₂h₂0)^2
    D₁₂ = (F₁h₁0 + F₂h₂0)*(F₃h₁0 + F₄h₂0)
    # D₂₂ = (F₃h₁0 + F₄h₂0)^2

    # `sqrtS` should be at `curvatureheight` because that is where `dFdθ` is evaluated
    T₁ = sqrtS*(1 + Rg[1,1])^2*(1 - R[2,2]*Rg[2,2])/(dFdθ*Rg[1,1]*D₁₁)
    # T₂ = sqrtS*(1 + Rg[2,2])^2*(1 - R[1,1]*Rg[1,1])/(dFdθ*Rg[2,2]*D₂₂)
    T₃ = sqrtS*(1 + Rg[1,1])*(1 + Rg[2,2])*R[2,1]/(dFdθ*D₁₂)
    T₄ = R[1,2]/R[2,1]

    # These are [^Pappert1983] terms divided by `-S`, the factor between Hy and Ez
    λv = -S₀*T₁
    λb = T₃*T₄
    λe = T₁

    return λv, λb, λe
end

@doc raw"""
    heightgains(z, ea₀, frequency, efconstants::ExcitationFactor; params=LMPParams())

Compute height gain functions at height `z` returned as the tuple `(fz, fy, fx)` where
eigenangle `ea₀` is referenced to the ground.

- `fz` is the height gain for the vertical electric field component ``Ez``.
- `fy` is the height gain for the transverse electric field component ``Ey``.
- `fx` is the height gain for the horizontal electric field component ``Ex``.
[^Pappert1983]

!!! note

    This function assumes that reflection coefficients are referenced to ``d = z = 0``.

See also: [`excitationfactorconstants`](@ref)

# References

[^Pappert1983]: R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low
    Frequency/Very Low Frequency) long path pulse program for antennas of arbitrary
    elevation and orientation,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891,
    Aug. 1983. Accessed: Jul. 04, 2018. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA133876.

[^Pappert1986]: R. A. Pappert and J. A. Ferguson, “VLF/LF mode conversion model calculations
    for air to air transmissions in the earth-ionosphere waveguide,” Radio Sci., vol. 21,
    no. 4, pp. 551–558, Jul. 1986, doi: 10.1029/RS021i004p00551.
"""
function heightgains(z, ea₀, frequency, efconstants::ExcitationFactor; params=LMPParams())
    C, C² = ea₀.cosθ, ea₀.cos²θ
    k = frequency.k
    @unpack F₁, F₂, F₃, F₄, Rg = efconstants
    @unpack earthradius, earthcurvature = params

    if earthcurvature
        # Precompute
        α = 2/earthradius
        expz = exp(z/earthradius)  # assumes reflection coefficients are referenced to `d = 0`

        qz = pow23(k/α)*(C² + α*z)  # (k/α)^(2/3)*(C² + α*z)

        h₁z, h₂z, dh₁z, dh₂z = modifiedhankel(qz)

        # Precompute
        F₁h₁z = F₁*h₁z
        F₂h₂z = F₂*h₂z

        # Height gain for Ez, also called f∥(z).
        fz = expz*(F₁h₁z + F₂h₂z)

        # Height gain for Ey, also called f⟂(z)
        fy = (F₃*h₁z + F₄*h₂z)

        # Height gain for Ex, also called g(z)
        # f₂ = 1/(1im*k) df₁/dz
        fx = expz/(1im*k*earthradius)*(F₁h₁z + F₂h₂z + earthradius*(F₁*dh₁z + F₂*dh₂z))
    else
        # BUG? I'm suspicious of the values, but would need to get modifiedhankel working
        # with big floats to check against earthradius ≈ Inf
        # Flat earth, [^Pappert1983] pg. 12--13
        expiz = cis(k*C*z)
        expmiz = cis(-k*C*z)
        fz = expiz + Rg[1,1]*expmiz
        fy = expiz + Rg[2,2]*expmiz
        fx = C*(expiz - Rg[1,1]*expmiz)
    end

    return fz, fy, fx
end

@doc raw"""
    modeterms(modeequation, tx::Emitter, rx::AbstractSampler; params=LMPParams())

Compute `tx` and `rx` height-gain and excitation factor products and `ExcitationFactor`
constants returned as the tuple `txterm, rxterm`.

The returned `txterm` is:
```math
λ_v \cos(γ) f_z(zₜ) + λ_b \sin(γ)\sin(ϕ) f_y(zₜ) + λ_e \sin(γ)\cos(ϕ) f_z(zₜ)
```
and `rxterm` is the height-gain function appropriate for `rx.fieldcomponent`:

| `fieldcomponent` |    f(zᵣ)       |
|:----------------:|:--------------:|
|      ``z``       |  ``-S₀⋅f_z``   |
|      ``y``       |  ``EyHy⋅f_y``  |
|      ``x``       |     ``-f_x``   |

# References

[^Pappert1976]: R. A. Pappert and L. R. Shockey, “Simplified VLF/LF mode conversion program
    with allowance for elevated, arbitrarily oriented electric dipole antennas,” Naval
    Electronics Laboratory Center, San Diego, CA, Interim Report 771, Oct. 1976. [Online].
    Available: http://archive.org/details/DTIC_ADA033412.

[^Morfitt1980]: D. G. Morfitt, “‘Simplified’ VLF/LF mode conversion computer programs:
    GRNDMC and ARBNMC,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-514, Jan. 1980.
    Accessed: Jan. 15, 2018. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA082695.
"""
function modeterms(modeequation, tx::Emitter, rx::AbstractSampler; params=LMPParams())
    @unpack ea, frequency, waveguide = modeequation
    @unpack ground = waveguide

    ea₀ = referencetoground(ea, params=params)
    S₀ = ea₀.sinθ

    frequency == tx.frequency ||
        throw(ArgumentError("`tx.frequency` and `modeequation.frequency` do not match"))

    zt = altitude(tx)
    zr = altitude(rx)
    rxfield = fieldcomponent(rx)

    # Transmit antenna orientation with respect to propagation direction
    # See [^Morfitt1980] pg. 22
    Sγ, Cγ = sincos(inclination(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(azimuth(tx))  # ϕ is measured from `x`

    t1 = Cγ
    t2 = Sγ*Sϕ
    t3 = Sγ*Cϕ

    dFdθ, R, Rg = solvedmodalequation(modeequation, params=params)
    efconstants = excitationfactorconstants(ea₀, R, Rg, frequency, ground, params=params)

    λv, λb, λe = excitationfactor(ea, dFdθ, R, efconstants, params=params)

    # Transmitter term
    fzt, fyt, fxt = heightgains(zt, ea₀, frequency, efconstants, params=params)
    txterm = λv*fzt*t1 + λb*fyt*t2 + λe*fxt*t3

    # Receiver term
    if zr == zt
        fzr, fyr, fxr = fzt, fyt, fxt
    else
        fzr, fyr, fxr = heightgains(zr, ea₀, frequency, efconstants, params=params)
    end

    # TODO: Handle multiple fields - maybe just always return all 3?
    if rxfield == Fields.Ez
        rxterm = -S₀*fzr
    elseif rxfield == Fields.Ey
        rxterm = efconstants.EyHy*fyr
    elseif rxfield == Fields.Ex
        rxterm = -fxr
    end

    return txterm, rxterm
end

"""
    modeterms(modeequation::ModeEquation, tx::Transmitter{VerticalDipole},
        rx::GroundSampler; params=LMPParams())

Specialized `modeterms` for the common case of `GroundSampler` and
`Transmitter{VerticalDipole}`.
"""
function modeterms(modeequation::ModeEquation, tx::Transmitter{VerticalDipole},
    rx::GroundSampler; params=LMPParams())

    @unpack ea, frequency, waveguide = modeequation
    @unpack ground = waveguide
    ea₀ = referencetoground(ea, params=params)
    S₀ = ea₀.sinθ

    frequency == tx.frequency ||
        throw(ArgumentError("`tx.frequency` and `modeequation.frequency` do not match"))

    rxfield = fieldcomponent(rx)

    dFdθ, R, Rg = solvedmodalequation(modeequation, params=params)
    efconstants = excitationfactorconstants(ea₀, R, Rg, frequency, ground, params=params)

    λv, λb, λe = excitationfactor(ea, dFdθ, R, efconstants, params=params)

    # Transmitter term
    # TODO: specialized heightgains for z = 0
    fz, fy, fx = heightgains(0.0, ea₀, frequency, efconstants, params=params)
    txterm = λv*fz

    # Receiver term
    if rxfield == Fields.Ez
        rxterm = -S₀*fz
    elseif rxfield == Fields.Ey
        rxterm = efconstants.EyHy*fy
    elseif rxfield == Fields.Ex
        rxterm = -fx
    end

    return txterm, rxterm
end

#==
Electric field calculation
==#

"""
    Efield(modes::Vector{EigenAngle}, waveguide::HomogeneousWaveguide, tx::Emitter,
           rx::AbstractSampler; params=LMPParams())

Compute the complex electric field by summing `modes` in `waveguide` with transmitter `tx`
and receiver `rx`.

See also: [`Efield!`](@ref)
"""
function Efield(modes::Vector{EigenAngle}, waveguide::HomogeneousWaveguide, tx::Emitter,
    rx::AbstractSampler; params=LMPParams())

    X = distance(rx, tx)
    Xlength = length(X)

    E = zeros(ComplexF64, Xlength)
    Efield!(E, modes, waveguide, tx, rx, params=params)

    return E
end

"""
    Efield!(E, modes::Vector{EigenAngle}, waveguide::HomogeneousWaveguide,
            tx::Emitter, rx::AbstractSampler; params=LMPParams())

Compute the complex electric field `E` in-place.

See also: [`Efield`](@ref)

# References

[^Morfitt1980]: D. G. Morfitt, “‘Simplified’ VLF/LF mode conversion computer programs:
    GRNDMC and ARBNMC,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-514, Jan. 1980.
    Accessed: Jan. 15, 2018. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA082695.

[^Pappert1983]: R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low
    Frequency/Very Low Frequency) long path pulse program for antennas of arbitrary
    elevation and orientation,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891,
    Aug. 1983. Accessed: Jul. 04, 2018. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA133876.
"""
function Efield!(E, modes::Vector{EigenAngle}, waveguide::HomogeneousWaveguide, tx::Emitter,
    rx::AbstractSampler; params=LMPParams())

    X = distance(rx, tx)
    length(E) == length(X) ||
        throw(ArgumentError("`E` must be same length as `rx` distances."))

    txpower = power(tx)
    frequency = tx.frequency
    k = frequency.k

    # Initialize `E` to zeros if necessary
    iszero(E) || fill!(E, 0)

    for ea in modes
        modeequation = PhysicalModeEquation(ea, frequency, waveguide)
        txterm, rxterm = modeterms(modeequation, tx, rx, params=params)

        S₀ = referencetoground(ea.sinθ, params=params)
        expterm = -k*(S₀ - 1)
        txrxterm = txterm*rxterm

        @inbounds for i in eachindex(E)
            E[i] += txrxterm*cis(expterm*X[i])
        end
    end

    Q = 0.6822408*sqrt(frequency.f*txpower)  # factor from lw_sum_modes.for
    # Q = Z₀/(4π)*sqrt(2π*txpower/10k)*k/2  # Ferguson and Morfitt 1981 eq (21), V/m, NOT uV/m!
    # Q *= 100 # for V/m to uV/m

    # TODO: Radiation resistance correction if zt > 0
    # See, e.g. Pappert Hitney 1989 TWIRE paper

    @inbounds for i in eachindex(E)
        E[i] *= Q/sqrt(abs(sin(X[i]/params.earthradius)))
    end

    # At transmitter (within 1 meter from it), E is complex NaN or Inf
    if X[1] < 1
        # Used in LWPC `lw_sum_modes.for`, but not sure where they got it
        # amplitude = 10log10(80*Q)
        E[1] = sqrt(80*Q) + 0.0im # == 10^(amplitude/20)
    end

    return E
end

"""
    Efield(modes::Vector{EigenAngle}, waveguide::HomogeneousWaveguide, tx::Emitter,
           rx::AbstractSampler{<:Real}; params=LMPParams())

Compute the scalar electric field when `AbstractSampler` has a (single) `Real` `distance`.
"""
function Efield(modes::Vector{EigenAngle}, waveguide::HomogeneousWaveguide, tx::Emitter,
    rx::AbstractSampler{<:Real}; params=LMPParams())

    # Unpack
    frequency = tx.frequency
    k = frequency.k

    txpower = power(tx)
    x = distance(rx, tx)

    Q = 0.6822408*sqrt(frequency.f*txpower)

    # At transmitter (within 1 meter from it), E is complex NaN or Inf
    if x < 1
        # Used in LWPC `lw_sum_modes.for`, but not sure where they got it
        # amplitude = 10log10(80*Q)
        E = sqrt(80*Q) + 0.0im # == 10^(amplitude/20)

        return E
    end

    E = zero(ComplexF64)
    for ea in modes
        modeequation = PhysicalModeEquation(ea, frequency, waveguide)
        txterm, rxterm = modeterms(modeequation, tx, rx, params=params)

        S₀ = referencetoground(ea.sinθ, params=params)
        expterm = -k*(S₀ - 1)
        txrxterm = txterm*rxterm

        E += txrxterm*cis(expterm*x)
    end

    E *= Q/sqrt(abs(sin(x/params.earthradius)))

    return E
end

"""
    Efield(waveguide::SegmentedWaveguide, wavefields_vec, adjwavefields_vec, tx::Emitter,
           rx::AbstractSampler; params=LMPParams())
"""
function Efield(waveguide::SegmentedWaveguide, wavefields_vec, adjwavefields_vec,
    tx::Emitter, rx::AbstractSampler; params=LMPParams())
    @unpack earthradius = params

    # Checks
    first(waveguide).distance == 0 ||
        throw(ArgumentError("The first `waveguide` segment should have `distance` 0.0."))
    length(waveguide) == length(wavefields_vec) == length(adjwavefields_vec) ||
        throw(ArgumentError("`wavefields_vec` and `adjwavefields_vec` must have the same"*
                            "length as `waveguide`."))

    X = distance(rx, tx)
    maxX = maximum(X)
    Xlength = length(X)
    E = Vector{ComplexF64}(undef, Xlength)

    frequency = tx.frequency
    k = frequency.k

    Q = 0.6822408*sqrt(frequency.f*tx.power)

    # Initialize
    J = length(waveguide)
    M = 0  # number of eigenangles in previous segment. Current segment is N
    xmtrfields = Vector{ComplexF64}(undef, 0)  # fields generated by transmitter
    previous_xmtrfields = similar(xmtrfields)  # fields saved from previous segment
    rcvrfields = similar(xmtrfields)  # fields at receiver location

    i = 1  # index of X
    for j = 1:J  # index of waveguide
        wvg = waveguide[j]
        wavefields = wavefields_vec[j]
        eas = eigenangles(wavefields)
        N = nummodes(wavefields)

        # Identify distance at beginning of segment
        segment_start = wvg.distance
        maxX < segment_start && break  # no farther X; break

        # Identify distance at end of segment
        if j < J
            segment_end = waveguide[j+1].distance
        else
            # last segment
            segment_end = typemax(typeof(segment_start))
        end

        # xmtrfields is for `Hy`
        resize!(xmtrfields, N)
        resize!(rcvrfields, N)
        if j > 1
            adjwavefields = adjwavefields_vec[j]
            prevwavefields = wavefields_vec[j-1]
            conversioncoeffs = modeconversion(prevwavefields, wavefields, adjwavefields,
                                              params=params)
        end

        # Calculate the mode terms (height gains and excitation factors) up to the current
        # segment
        for n = 1:N
            modeequation = PhysicalModeEquation(eas[n], frequency, wvg)
            txterm, rxterm = modeterms(modeequation, tx, rx, params=params)
            if j == 1
                # Transmitter exists only in the transmitter slab (obviously)
                xmtrfields[n] = txterm
            else
                # Otherwise, mode conversion of transmitted fields
                xmtrfields_sum = zero(eltype(xmtrfields))
                for m = 1:M
                    xmtrfields_sum += previous_xmtrfields[m]*conversioncoeffs[m,n]
                end
                xmtrfields[n] = xmtrfields_sum
            end

            rcvrfields[n] = xmtrfields[n]*rxterm
        end

        # Calculate E at each distance in the current waveguide segment
        while X[i] < segment_end
            x = X[i] - segment_start
            factor = Q/sqrt(abs(sin(X[i]/earthradius)))

            totalfield = zero(eltype(E))
            for n = 1:N
                S₀ = referencetoground(eas[n].sinθ, params=params)
                totalfield += rcvrfields[n]*cis(-k*x*(S₀ - 1))*factor
            end

            E[i] = totalfield
            i += 1
            i > Xlength && break
        end

        # If we've reached the end of the current segment and there are more segments,
        # prepare for next segment
        if j < J
            # End of current slab
            x = segment_end - segment_start

            resize!(previous_xmtrfields, N)
            for n = 1:N
                S₀ = referencetoground(eas[n].sinθ, params=params)

                # Excitation factors at end of slab
                xmtrfields[n] *= cis(-k*x*(S₀ - 1))
                previous_xmtrfields[n] = xmtrfields[n]
            end
            M = N  # set previous number of modes
        end
    end

    # At transmitter (within 1 meter from it), E is complex NaN or Inf
    if X[1] < 1
        # Used in LWPC `lw_sum_modes.for`, but not sure where they got it
        # amplitude = 10log10(80*Q)
        E[1] = sqrt(80*Q) + 0.0im # == 10^(amplitude/20)
    end

    return E
end

########

"""
radiation resistance correction factor for when zt isn't 0.

From lw_sum_modes.for
but could also see Pappert and Hitney 1989 TWIRE paper
"""
function radiationresistance(k, Cγ, zt)
    x = 2*k*zt
    sinx, cosx = sincos(x)
    xt1 = 3*(sinx - x*cosx)/x^3
    xt2 = (xt1 - 3*sinx/x)/2
    xt3 = sqrt(2/(1 + xt2 + (xt1 - xt2)*Cγ^2))

    return xt3
end

#==
TODO: LF corrections/fixes

Pappert 1981 LF Daytime Earth Ionosphere...
finds that the linear combination of modified Hankel functions of order one third used to
represent the height gain at the ground is in fact incorrect at the ground for modes which
are highly earth detached (pg 7-8). To correct for this, they just throw away earth
effects altogether when the condition
Re(2im*(k/α)*(C²ₕ - α*H)^(3/2)/3)) > 12.4
is met. The value of 12.4 was found by trial and error and requires the degree of evanescence
at the ground to be of the order of several times 10⁻⁶. When this condition is met, the
plane wave reflection coefficients (for ground) become eq (2) and (3).
There are additional equation replacements on page 10.
Otherwise, this is following the math from Pappert & Shockey 71. This paper (81) explicitly
identifies where some of the angles should be referenced.
==#

#==
TODO: ELF corrections/fixes

Pappert Shockey 1971 WKB Mode Summing...
changes height gain functions for ELF, pg. 9

Switches to flat earth at imaginary angles less than -10° (see LWPC or Pappert 1983 ELF-VLF)

Pappert Shockey 1971 pg 9 or Pappert 198e pg 12
elf_heightgains()
==#
