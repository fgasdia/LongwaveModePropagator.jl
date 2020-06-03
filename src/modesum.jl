@with_kw struct ExcitationFactor{T} @deftype T
    F₁
    F₂
    F₃
    F₄
    h₁0
    h₂0
end

################
# Excitation and Mode Sum
################

"""
    excitationfactorconstants(ea, R, Rg, frequency, ground)

Return an `ExcitationFactor` struct used in calculating height gains.

Based on Morfitt 1980, Pappert Shockey 1971, and Pappert Shockey 1976 (this last one has H=0)

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Ferguson1981] J. A. Ferguson and D. G. Morfitt, “WKB mode summing program for dipole antennas of arbitrary orientation and elevation for VLF/LF propagation,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-697, Oct. 1981.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function excitationfactorconstants(ea::EigenAngle, R, Rg, frequency::Frequency, ground::Ground)
    S², C² = ea.sin²θ, ea.cos²θ
    k, ω = frequency.k, frequency.ω

    # `ea` is at height `CURVATURE_HEIGHT`. See, e.g. Pappert1971 pg 8

    # Precompute
    α = 2/Rₑ
    αH = α*CURVATURE_HEIGHT
    koα = k/α
    cbrtkoα = cbrt(koα)
    koα23 = cbrtkoα^2  # (k/α)^(2/3)
    αok23 = inv(koα23)  # (α/k)^(2/3)

    q₀ = koα23*(C² - αH)  # (k/α)^(2/3)*(C² - αH)

    # XXX: `modifiedhankel` dominates runtime of this function
    h₁0, h₂0, h₁p0, h₂p0 = modifiedhankel(q₀)

    H₁0 = h₁p0 + αok23*h₁0/2
    H₂0 = h₂p0 + αok23*h₂0/2

    n₀² = 1 - αH  # modified index of refraction (free space) squared
    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*ϵ₀))  # ground index of refraction

    # Precompute
    n₀²oNg² = n₀²/Ng²
    sqrtNg²mS² = sqrt(Ng² - S²)

    cbrtkoαsqrtNg²mS²h₂0 = cbrtkoα*sqrtNg²mS²*h₂0
    cbrtkoαsqrtNg²mS²h₁0 = cbrtkoα*sqrtNg²mS²*h₁0

    F₁ = -H₂0 + im*n₀²oNg²*cbrtkoαsqrtNg²mS²h₂0
    F₂ = H₁0 - im*n₀²oNg²*cbrtkoαsqrtNg²mS²h₁0
    F₃ = -h₂p0 + im*cbrtkoαsqrtNg²mS²h₂0
    F₄ = h₁p0 - im*cbrtkoαsqrtNg²mS²h₁0

    # ey/hy; polarization ratio; Normalizes y component of H to unity at thr ground.
    # Sometimes called `f` in papers
    # f0fr = T₂/(T₃*T₄) = T₃/T₁
    # LWPC uses the following criteria for choosing
    # if abs2(1 - R[1,1]*Rg[1,1]) > abs2(1 - R[2,2]*Rg[2,2])
    #     f0fr = (1 + Rg[2,2])*(1 - R[1,1]*Rg[1,1])/((1 + Rg[1,1])*R[1,2]*Rg[2,2])
    # else
    #     f0fr = (1 + Rg[2,2])*R[2,1]*Rg[1,1]/((1 + Rg[1,1])*(1 - R[2,2]*Rg[2,2]))
    # end

    return ExcitationFactor(F₁, F₂, F₃, F₄, h₁0, h₂0)
end

"""
    heightgains(z, ea, frequency, efconstants)

Calculate heightgains at `z`.

This function assumes the reference height for the reflection coefficients is ``d=0``.

See also: [`excitationfactor`](@ref), [`excitationfactorconstants`](@ref)

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function heightgains(z, ea::EigenAngle, frequency::Frequency, efconstants::ExcitationFactor)
    C² = ea.cos²θ
    k = frequency.k
    @unpack F₁, F₂, F₃, F₄, h₁0, h₂0 = efconstants

    # Precompute
    α = 2/Rₑ
    koα23 = pow23(k/α)  # (k/α)^(2/3)
    expzoRₑ = exp(z/Rₑ)  # assumes `d = 0`

    qz = koα23*(C² - α*(CURVATURE_HEIGHT - z))

    # XXX: `modifiedhankel` dominates this functions runtime
    h₁z, h₂z, h₁pz, h₂pz = modifiedhankel(qz)

    # Precompute
    F₁h₁0 = F₁*h₁0
    F₂h₂0 = F₂*h₂0
    F₁h₁z = F₁*h₁z
    F₂h₂z = F₂*h₂z

    # From Ferguson1981 pg 10 or Pappert1971 pg 6, 8:
    # Height gain for Ez, also called f∥(z)
    fz = expzoRₑ*(F₁h₁z + F₂h₂z)

    # Height gain for Ex, also called g(z)
    # f₂ = 1/(im*k) df₁/dz
    fx = -im*expzoRₑ/(Rₑ*k)*(F₁h₁z + F₂h₂z + Rₑ*(F₁*h₁pz + F₂*h₂pz))

    # Height gain for Ey, also called f⟂(z)
    fy = (F₃*h₁z + F₄*h₂z)

    return fz, fx, fy
end

"""
    excitationfactor(ea, dFdθ, R, Rg, component)

Calculate the excitation factor for electric field `component`.

These excitation factors are used in conjunction with the function [`heightgains`](@ref).

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Ferguson1981] J. A. Ferguson and D. G. Morfitt, “WKB mode summing program for dipole antennas of arbitrary orientation and elevation for VLF/LF propagation,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-697, Oct. 1981.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center, San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function excitationfactor(ea::EigenAngle, dFdθ, R, Rg, efconstants::ExcitationFactor, component::FieldComponent)
    S, S², C² = ea.sinθ, ea.sin²θ, ea.cos²θ
    sqrtS = sqrt(S)

    @unpack F₁, F₂, F₃, F₄, h₁0, h₂0 = efconstants

    # TODO: Beter incorporation with other functions?
    # Calculate height gains for `z=d=0`
    fz = (F₁*h₁0 + F₂*h₂0)
    fy = (F₃*h₁0 + F₄*h₂0)

    D₁₁ = fz^2
    D₁₂ = fz*fy
    D₂₂ = fy^2

    # `S` is at `CURVATURE_HEIGHT`, specified in e.g.
    # D. G. Morfitt, ``'Simplified' VLF/LF mode conversion...,'' NOSC/TR-514, 1980, pg 19.
    T₁ = sqrtS*(1 + Rg[1,1])^2*(1 - R[2,2]*Rg[2,2])/(dFdθ*Rg[1,1]*D₁₁)
    T₂ = sqrtS*(1 + Rg[2,2])^2*(1 - R[1,1]*Rg[1,1])/(dFdθ*Rg[2,2]*D₂₂)
    T₃ = sqrtS*(1 + Rg[1,1])*(1 + Rg[2,2])*R[2,1]/(dFdθ*D₁₂)
    T₄ = R[1,2]/R[2,1]

    ST₁ = S*T₁
    ST₃ = S*T₃
    if component == FC_Ez
        λv = S²*T₁
        λb = -ST₃*T₄
        λe = -ST₁
    elseif component == FC_Ey
        λv = -ST₃
        λb = T₂
        λe = T₃
    elseif component == FC_Ex
        λv = ST₁
        λb = -T₃*T₄
        λe = -T₁
    end

    return λv, λb, λe
end

function modeterms(ea, frequency, waveguide, emitter_orientation, sampler_orientation)
    # Unpack
    Sγ, Cγ, Sϕ, Cϕ, zt = emitter_orientation
    rxcomponent, zr = sampler_orientation

    dFdθ, R, Rg = solvemodalequationdθ(ea, frequency, waveguide)
    efconstants = excitationfactorconstants(ea, R, Rg, frequency, waveguide.ground)

    # fz0, fx0, fy0 = heightgains(0, ea, frequency, efconstants)
    λv, λb, λe = excitationfactor(ea, dFdθ, R, Rg, efconstants, rxcomponent)

    # Transmitter term
    fz_t, fx_t, fy_t = heightgains(zt, ea, frequency, efconstants)
    # λv, λe, λb = excitationfactor(ea, dFdθ, R, Rg, efconstants, rxcomponent)
    xmtrterm = λv*fz_t*Cγ + λb*fy_t*Sγ*Sϕ + λe*fx_t*Sγ*Cϕ

    # Receiver term
    # TODO: Handle multiple components
    fz_r, fx_r, fy_r = heightgains(zr, ea, frequency, efconstants)
    if rxcomponent == FC_Ez
        rcvrterm = fz_r
    elseif rxcomponent == FC_Ex
        rcvrterm = fx_r
    elseif rxcomponent == FC_Ey
        rcvrterm = fy_r
    end

    return xmtrterm, rcvrterm
end

"""
    Efield!(E, X, modes, waveguide, tx, rx)

Calculate the complex electric field, amplitude, and phase at a distance `x` from transmitter `tx`.

`modes` is a collection of `EigenAngles` for the earth-ionosphere waveguide with the parameters
`modeparams`. Emitter `tx` specifies the transmitting antenna position, orientation, and
radiated power, and `rx` specifies the field component of interest.

NOTE: this returns modesum without `x`. To get correct values, need to raise sum to x power

# References


"""
function Efield!(E::AbstractVector{<:Complex}, modes, waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler)
    X = distance(rx,tx)
    @assert length(E) == length(X)  # or boundscheck

    # TODO: special function for vertical component, transmitter, and at ground
    # TODO: Special Efield() for point measurement

    frequency = tx.frequency

    txpower = power(tx)
    zt = altitude(tx)
    k = frequency.k
    zr = altitude(rx)
    rxcomponent = fieldcomponent(rx)

    # Transmit dipole antenna orientation with respect to propagation direction
    # See Morfitt 1980 pg 22
    # TODO: Confirm alignment of coord system and magnetic field
    Sγ, Cγ = sincos(π/2 - elevation(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(azimuth(tx))  # ϕ is measured from `x`

    emitter_orientation = (Sγ=Sγ, Cγ=Cγ, Sϕ=Sϕ, Cϕ=Cϕ, zt=zt)
    sampler_orientation = (rxcomponent=rxcomponent, zr=zr)

    # Initialize E if necessary
    iszero(E) || fill!(E, 0)

    for ea in modes
        xmtrterm, rcvrterm = modeterms(ea, frequency, waveguide, emitter_orientation, sampler_orientation)

        # Precalculate
        S₀ = referencetoground(ea.sinθ)
        expterm = -k*(S₀ - 1)
        xmtr_rcvr_term = xmtrterm*rcvrterm

        @inbounds for i in eachindex(E)
            E[i] += xmtr_rcvr_term*cis(expterm*X[i])
        end
    end

    Q = 682.2408*sqrt(frequency.f/1000*txpower/1000)  # in lw_sum_modes.for
    # Q = Z₀/(4π)*sqrt(2π*txpower/10k)*k/2  # Ferguson and Morfitt 1981 eq (21), V/m, NOT uV/m!
    # Q *= 100 # for V/m to uV/m

    # TODO: Radiation resistance correction if zt > 0

    @inbounds for i in eachindex(E)
        E[i] *= Q/sqrt(abs(sin(X[i]/Rₑ)))
    end

    return nothing
end

function Efield(
    modes,
    waveguide::HomogeneousWaveguide,
    tx::Emitter,
    rx::AbstractSampler
    )

    frequency = tx.frequency

    X = distance(rx, tx)
    Xlength = length(X)

    Etype = eltype(eltype(modes))
    E = zeros(Etype, Xlength)
    phase = Array{real(Etype)}(undef, Xlength)
    amp = Array{real(Etype)}(undef, Xlength)

    Efield!(E, modes, waveguide, tx, rx)

    @inbounds for i in eachindex(E)
        e = E[i]
        if
        phase[i] = angle(e)  # ranges between -π:π rad
        amp[i] = 10log10(abs2(e))  # == 20log10(abs(E))
    end

    unwrap!(phase)

    return E, phase, amp
end


function Efield!(E::AbstractVector{<:Complex}, waveguide::SegmentedWaveguide, tx::Emitter, rx::AbstractSampler)
    # catch if `waveguide` only has a single segment
    num_segments = length(waveguide)
    num_segments == 1 && return Efield(modes, only(waveguide), tx, rx)

    X = distance(rx, tx)
    @assert length(E) == length(X)  # or boundscheck

    # TODO: special function for vertical component, transmitter, and at ground
    # TODO: Special Efield() for point measurement

    frequency = tx.frequency

    txpower = power(tx)
    zt = altitude(tx)
    k = frequency.k
    zr = altitude(rx)
    rxcomponent = fieldcomponent(rx)

    zs = range(TOPHEIGHT, zero(TOPHEIGHT), length=257)

    # Transmit dipole antenna orientation with respect to propagation direction
    # See Morfitt 1980 pg 22
    # TODO: Confirm alignment of coord system and magnetic field
    Sγ, Cγ = sincos(π/2 - elevation(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(azimuth(tx))  # ϕ is measured from `x`

    emitter_orientation = (Sγ=Sγ, Cγ=Cγ, Sϕ=Sϕ, Cϕ=Cϕ, zt=zt)
    sampler_orientation = (rxcomponent=rxcomponent, zr=zr)

    # Initialize E if necessary
    iszero(E) || fill!(E, 0)



    origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
    origcoords .= deg2rad.(origcoords)
    tolerance = 1e-8


    # Transmitter segment is special
    # TODO: check transmitter slab is the first slab
    txwaveguide = first(waveguide)

    modes = findmodes(origcoords, frequency, txwaveguide, tolerance)

    wavefields = Wavefields(modes, zs)
    adjwavefields = Wavefields(modes, zs)

    calculate_wavefields!(wavefields, adjwavefields, frequency, waveguide)
    prevwavefields = copy(wavefields)

    # S₀t = referencetoground.(getindex.(modes, :sinθ))

    segment_range = waveguide[2].distance
    xi = 1
    for ea in modes
        xmtrterm, rcvrterm = modeterms(ea, frequency, txwaveguide, emitter_orientation, sampler_orientation)

        # Precalculate
        S₀ = referencetoground(ea.sinθ)
        expterm = -k*(S₀ - 1)
        xmtr_rcvr_term = xmtrterm*rcvrterm

        while X[xi] <= segment_range
            E[xi] += xmtr_rcvr_term*cis(expterm*X[xi])
            xi += 1
        end
    end


    lastj = lastindex(waveguide)
    for j in 2:num_segments
        segment_range = j == lastj ? oftype(waveguide[j].distance, 40000) : waveguide[j+1].distance

        modes = findmodes(origcoords, frequency, waveguide[j], tolerance)

        wavefields = Wavefields(modes, zs)
        adjwavefields = Wavefields(modes, zs)

        calculate_wavefields!(wavefields, adjwavefields, frequency, waveguide)

        a = modeconversion(prevwavefields, wavefields, adjwavefields)


        for ea in modes

            xmtrterm, rcvrterm = modeterms(ea, frequency, txwaveguide, emitter_orientation, sampler_orientation)

            # Precalculate
            S₀ = referencetoground(ea.sinθ)
            expterm = -k*(S₀ - 1)

            if rxcomponent != FC_Ez
                rcvrterm *= S₀t/S₀
            end

            while X[xi] <= segment_range
                @inbounds E[xi] += xmtr_rcvr_term*cis(expterm*X[xi])
                xi += 1
            end

            prevwavefields, wavefields = wavefields, prevwavefields
        end
        xi == lastindex(X) && break
    end

    Q = 682.2408*sqrt(frequency.f/1000*txpower/1000)  # in lw_sum_modes.for
    # Q = Z₀/(4π)*sqrt(2π*txpower/10k)*k/2  # Ferguson and Morfitt 1981 eq (21), V/m, NOT uV/m!
    # Q *= 100 # for V/m to uV/m

    # TODO: Radiation resistance correction if zt > 0

    @inbounds for i in eachindex(E)
        E[i] *= Q/sqrt(abs(sin(X[i]/Rₑ)))
    end

    return nothing
end






#==
"""
Pappert Shockey 1976
"""
function Efield(modes, waveguide::SegmentedWaveguide, tx::Emitter, rx::GroundSampler)




    # Morfitt 1980
    # `k` is mode index -> now `i`
    # `n` is field index 1, 2, 3
    # so `j` must be a slab index

==#

########
# Utility functions

function unwrap!(phasearray::AbstractVector)
    @inbounds for i in 2:length(phasearray)
        d = phasearray[i] - phasearray[i-1]
        if d >= π
            d -= 2π
        elseif d < -π
            d += 2π
        end
        phasearray[i] = phasearray[i-1] + d
    end
    return nothing
end

# see, e.g. PS71 pg 11
function referencetoground(ea::EigenAngle)
    return EigenAngle(asin(ea.sinθ/sqrt(1 - 2/Rₑ*CURVATURE_HEIGHT)))
end
referencetoground(x::Number) = x/sqrt(1 - 2/Rₑ*CURVATURE_HEIGHT)

"""
    pow23(x)

Calculate `x^(2/3)` relatively efficiently for `Real` or `Complex` `x`.
"""
pow23(x::Real) = cbrt(x)^2
pow23(z::Complex) = exp(2/3*log(z))

########

"""
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
"""

"""
TODO: ELF corrections/fixes

Pappert Shockey 1971 WKB Mode Summing...
changes height gain functions for ELF, pg. 9

Switches to flat earth at imaginary angles less than -10° (see LWPC or Pappert 1983 ELF-VLF)
"""

"""
Pappert Shockey 1971 pg 9 or Pappert 198e pg 12
"""
function elf_heightgains()

end
