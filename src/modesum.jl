@with_kw struct ExcitationFactor{T} @deftype T
    F₁
    F₂
    F₃
    F₄
    h₁0
    h₂0
    f0fr
end

################
# Excitation and Mode Sum
################

"""
    excitationfactorconstants(ea, R, Rg, frequency, ground)

Return an `ExcitationFactor` struct used in calculating height gains.

!!! note

    Eigenangle `ea` should be referenced to `CURVATURE_HEIGHT`. See, e.g.
    [^Pappert1971], page 8.

Based on Morfitt 1980, Pappert Shockey 1971, and Pappert Shockey 1976 (this last one has H=0)

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for
    VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics
    Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Ferguson1981] J. A. Ferguson and D. G. Morfitt, “WKB mode summing program for
    dipole antennas of arbitrary orientation and elevation for VLF/LF propagation,”
    Naval Ocean Systems Center, San Diego, CA, NOSC/TR-697, Oct. 1981.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF
    (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for
    Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center,
    San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function excitationfactorconstants(ea, R, Rg, frequency, ground)
    S², C² = ea.sin²θ, ea.cos²θ
    k, ω = frequency.k, frequency.ω

    # This function could be further optimized by reducing repeated computations,
    # but <100 ns can be gained and this function is not in a hot loop.
    # Clarity outways efficiency here.

    # Precompute
    α = 2/EARTHRADIUS
    αH = α*CURVATURE_HEIGHT
    t0 = (α/k)^(2/3)/2

    q₀ = (k/α)^(2/3)*(C² - αH)
    h₁0, h₂0, h₁p0, h₂p0 = modifiedhankel(q₀)

    H₁0 = h₁p0 + t0*h₁0
    H₂0 = h₂p0 + t0*h₂0

    n₀² = 1 - αH  # modified index of refraction (free space) squared
    Ng² = complex(ground.ϵᵣ, -ground.σ/(ω*E0))  # ground index of refraction

    # Precompute
    t1 = 1im*cbrt(k/α)*sqrt(Ng² - S²)

    F₁ = -H₂0 + (n₀²/Ng²)*t1*h₂0
    F₂ = H₁0 - (n₀²/Ng²)*t1*h₁0
    F₃ = -h₂p0 + t1*h₂0
    F₄ = h₁p0 - t1*h₁0

    # ey/hy; polarization ratio; Normalizes y component of H to unity at thr ground.
    # Sometimes called `f` in papers
    # f0fr = T₂/(T₃*T₄) = T₃/T₁
    # LWPC uses the following criteria for choosing
    if abs2(1 - R[1,1]*Rg[1,1]) > abs2(1 - R[2,2]*Rg[2,2])
        f0fr = (1 + Rg[2,2])*(1 - R[1,1]*Rg[1,1])/((1 + Rg[1,1])*R[1,2]*Rg[2,2])
    else
        f0fr = (1 + Rg[2,2])*R[2,1]*Rg[1,1]/((1 + Rg[1,1])*(1 - R[2,2]*Rg[2,2]))
    end

    return ExcitationFactor(F₁, F₂, F₃, F₄, h₁0, h₂0, f0fr)
end

"""
    heightgains(z, ea, frequency, efconstants)

Calculate heightgains at `z`.

This function assumes the reference height for the reflection coefficients is ``d=0``.

See also: [`excitationfactor`](@ref), [`excitationfactorconstants`](@ref)

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for
    VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics
    Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF
    (Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for
    Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center,
    San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function heightgains(z, ea, frequency, efconstants::ExcitationFactor)
    C² = ea.cos²θ
    k = frequency.k
    @unpack F₁, F₂, F₃, F₄, h₁0, h₂0 = efconstants

    # Precompute
    α = 2/EARTHRADIUS
    t1 = exp(z/EARTHRADIUS)  # assumes `d = 0`

    qz = (k/α)^(2/3)*(C² - α*(CURVATURE_HEIGHT - z))

    h₁z, h₂z, h₁pz, h₂pz = modifiedhankel(qz)

    # Precompute
    F₁h₁0 = F₁*h₁0
    F₂h₂0 = F₂*h₂0
    F₁h₁z = F₁*h₁z
    F₂h₂z = F₂*h₂z

    # From Ferguson1981 pg 10 or Pappert1971 pg 6, 8:
    # Height gain for Ez, also called f∥(z)
    fz = t1*(F₁h₁z + F₂h₂z)

    # Height gain for Ex, also called g(z)
    # f₂ = 1/(im*k) df₁/dz
    fx = -1im*t1/(EARTHRADIUS*k)*(F₁h₁z + F₂h₂z + EARTHRADIUS*(F₁*h₁pz + F₂*h₂pz))

    # Height gain for Ey, also called f⟂(z)
    fy = (F₃*h₁z + F₄*h₂z)

    return fz, fx, fy
end

"""
    excitationfactor(ea, dFdθ, R, Rg, field)

Calculate the excitation factor for electric `field`.

These excitation factors are used in conjunction with the function
[`heightgains`](@ref).

!!! note

    Eigen angle `ea` should be referenced to `CURVATURE_HEIGHT`. This function
    internally references `ea` to ground where necessary.

# References

[^Pappert1971] R. A. Pappert and L. R. Shockey, “WKB Mode Summing Program for
VLF/ELF Antennas of Arbitrary Length, Shape and Elevation,” Naval Electronics
Lab Center, San Diego, CA, NELC-IR-713, M402, Jun. 1971.

[^Ferguson1981] J. A. Ferguson and D. G. Morfitt, “WKB mode summing program for
dipole antennas of arbitrary orientation and elevation for VLF/LF propagation,”
Naval Ocean Systems Center, San Diego, CA, NOSC/TR-697, Oct. 1981.

[^Morfitt1980] D. G. Morfitt, “‘Simplified’ VLF/LF mode conversion computer
programs: GRNDMC and ARBNMC,” Naval Ocean Systems Center, San Diego, CA,
NOSC/TR-514, Jan. 1980.

[^Pappert1983] R. A. Pappert, L. R. Hitney, and J. A. Ferguson, “ELF/VLF
(Extremely Low Frequency/Very Low Frequency) Long Path Pulse Program for
Antennas of Arbitrary Elevation and Orientation.,” Naval Ocean Systems Center,
San Diego, CA, NOSC/TR-891, Aug. 1983.
"""
function excitationfactor(ea, dFdθ, R, Rg, efconstants::ExcitationFactor,
    field::Fields.Field)

    S, S², C² = ea.sinθ, ea.sin²θ, ea.cos²θ
    sqrtS = sqrt(S)

    S₀ = referencetoground(S)

    @unpack F₁, F₂, F₃, F₄, h₁0, h₂0 = efconstants

    # Calculate height gains for `z = d = 0`
    fz = (F₁*h₁0 + F₂*h₂0)
    fy = (F₃*h₁0 + F₄*h₂0)

    D₁₁ = fz^2
    D₁₂ = fz*fy
    D₂₂ = fy^2

    # `S` is at `CURVATURE_HEIGHT`, specified in e.g. [^Morfitt1980]
    T₁ = sqrtS*(1 + Rg[1,1])^2*(1 - R[2,2]*Rg[2,2])/(dFdθ*Rg[1,1]*D₁₁)
    T₂ = sqrtS*(1 + Rg[2,2])^2*(1 - R[1,1]*Rg[1,1])/(dFdθ*Rg[2,2]*D₂₂)
    T₃ = sqrtS*(1 + Rg[1,1])*(1 + Rg[2,2])*R[2,1]/(dFdθ*D₁₂)
    T₄ = R[1,2]/R[2,1]

    # Here eigenangle is at the ground (as is everything else)
    S₀T₁ = S₀*T₁
    S₀T₃ = S₀*T₃
    if field == Fields.Ez
        λv = S₀*S₀T₁
        λb = -S₀T₃*T₄
        λe = -S₀T₁
    elseif field == Fields.Ey
        λv = -S₀T₃
        λb = T₂
        λe = T₃
    elseif field == Fields.Ex
        λv = S₀T₁
        λb = -T₃*T₄
        λe = -T₁
    end

    return λv, λb, λe
end

"""
    modeterms(ea, frequency, waveguide, tx, rx)

!!! note

    Eigenangle `ea` should be referenced to `CURVATURE_HEIGHT`.
"""
function modeterms(modeequation, tx::Emitter, rx::AbstractSampler)
    @unpack ea, frequency, waveguide = modeequation

    zt = altitude(tx)
    zr = altitude(rx)
    rxfield = fieldcomponent(rx)

    # Transmit antenna orientation with respect to propagation direction
    # See Morfitt 1980 pg 22
    # TODO: Confirm alignment of coord system and magnetic field
    Sγ, Cγ = sincos(inclination(tx))  # γ is measured from vertical
    Sϕ, Cϕ = sincos(azimuth(tx))  # ϕ is measured from `x`

    t1 = Cγ
    t2 = Sγ*Sϕ
    t3 = Sγ*Cϕ

    dFdθ, R, Rg = solvemodalequation(ea, modeequation, Dθ())
    efconstants = excitationfactorconstants(ea, R, Rg, frequency, waveguide.ground)

    λv, λb, λe = excitationfactor(ea, dFdθ, R, Rg, efconstants, rxfield)

    # Transmitter term
    fz_t, fx_t, fy_t = heightgains(zt, ea, frequency, efconstants)
    xmtrterm = λv*fz_t*t1 + λb*fy_t*t2 + λe*fx_t*t3

    # Receiver term
    if zr == zt
        fz_r, fx_r, fy_r = fz_t, fx_t, fy_t
    else
        fz_r, fx_r, fy_r = heightgains(zr, ea, frequency, efconstants)
    end

    # TODO: Handle multiple fields - maybe just always return all 3?
    if rxfield == Fields.Ez
        rcvrterm = fz_r
    elseif rxfield == Fields.Ex
        rcvrterm = fx_r
    elseif rxfield == Fields.Ey
        rcvrterm = fy_r
    end

    return xmtrterm, rcvrterm, efconstants
end

"""
    modeterms()

Specialized `modeterms` for the common case of `GroundSampler` and
`Transmitter{VerticalDipole}`.
"""
function modeterms(modeequation::ModeEquation, tx::Transmitter{VerticalDipole},
    rx::GroundSampler)

    @unpack ea, frequency, waveguide = modeequation

    rxfield = fieldcomponent(rx)

    dFdθ, R, Rg = solvemodalequation(modeequation, Dθ())
    efconstants = excitationfactorconstants(ea, R, Rg, frequency, waveguide.ground)

    λv, λb, λe = excitationfactor(ea, dFdθ, R, Rg, efconstants, rxfield)

    # Transmitter term
    fz, fx, fy = heightgains(0.0, ea, frequency, efconstants)
    xmtrterm = λv*fz

    # Receiver term
    # TODO: Handle multiple fields - maybe just always return all 3?
    if rxfield == Fields.Ez
        rcvrterm = fz
    elseif rxfield == Fields.Ex
        rcvrterm = fx
    elseif rxfield == Fields.Ey
        rcvrterm = fy
    end

    return xmtrterm, rcvrterm, efconstants
end

###

"""
    Efield(modes, waveguide, tx, rx)

Calculate the complex electric field by summing `modes` in `waveguide` with
transmitter `tx` and receiver `rx`.
"""
function Efield(
    modes::Vector{EigenAngle},
    waveguide::HomogeneousWaveguide,
    tx::Emitter,
    rx::AbstractSampler
    )

    X = distance(rx, tx)
    Xlength = length(X)

    E = zeros(ComplexF64, Xlength)
    Efield!(E, X, modes, waveguide, tx, rx)

    return E
end

"""
    Efield!(E, X, modes, waveguide, tx, rx)

Calculate the complex electric field at a distance `x` from transmitter `tx`.

`modes` is a collection of `EigenAngles` for the earth-ionosphere waveguide with the parameters
`modeparams`. Emitter `tx` specifies the transmitting antenna position, orientation, and
radiated power, and `rx` specifies the field component of interest.
"""
function Efield!(
    E,
    X::AbstractVector,
    modes::AbstractVector,
    waveguide::HomogeneousWaveguide,
    tx::Emitter,
    rx::AbstractSampler
    )
    @assert length(E) == length(X) "`E` and `X` must be same length"

    frequency = tx.frequency

    txpower = power(tx)
    k = frequency.k

    # Initialize E if necessary
    iszero(E) || fill!(E, 0)

    for ea in modes
        modeequation = PhysicalModeEquation(ea, frequency, waveguide)
        xmtrterm, rcvrterm = modeterms(modeequation, tx, rx)

        # Precalculate
        S₀ = referencetoground(ea.sinθ)
        expterm = -k*(S₀ - 1)
        xmtr_rcvr_term = xmtrterm*rcvrterm

        @inbounds for i in eachindex(E)
            E[i] += xmtr_rcvr_term*cis(expterm*X[i])
        end
    end

    Q = 0.6822408*sqrt(frequency.f*txpower)  # factor from lw_sum_modes.for
    # Q = Z₀/(4π)*sqrt(2π*txpower/10k)*k/2  # Ferguson and Morfitt 1981 eq (21), V/m, NOT uV/m!
    # Q *= 100 # for V/m to uV/m

    # TODO: Radiation resistance correction if zt > 0
    # See, e.g. Pappert Hitney 1989 TWIRE paper

    @inbounds for i in eachindex(E)
        E[i] *= Q/sqrt(abs(sin(X[i]/EARTHRADIUS)))
    end

    return nothing
end

"""
    Efield(modes, waveguide, tx, rx)

Return the single scalar electric field when `AbstractSampler` has a (single)
`Real` distance.
"""
function Efield(
    modes::Vector{EigenAngle},
    waveguide::HomogeneousWaveguide,
    tx::Emitter,
    rx::AbstractSampler{R}
    ) where {R<:Real}

    # Unpack
    frequency = tx.frequency
    k = frequency.k

    txpower = power(tx)
    x = distance(rx, tx)

    E = zero(T)
    for ea in modes
        modeequation = PhysicalModeEquation(ea, frequency, waveguide)
        xmtrterm, rcvrterm = modeterms(modeequation, tx, rx)

        S₀ = referencetoground(ea.sinθ)
        expterm = -k*(S₀ - 1)
        xmtr_rcvr_term = xmtrterm*rcvrterm

        E += xmtr_rcvr_term*cis(expterm*x)
    end

    return E
end

function Efield(
    waveguide,
    wavefields_vec::Vector{Wavefields{T}},
    adjwavefields_vec::Vector{Wavefields{T}},
    tx,
    rx
    ) where {T}

    X = distance(rx, tx)
    maxX = maximum(X)
    Xlength = length(X)
    E = Vector{T}(undef, Xlength)

    frequency = tx.frequency
    k = frequency.k

    Q = 0.6822408*sqrt(frequency.f*tx.power)

    # Initialize
    numsegments = length(waveguide)
    previous_eacount = 0  # number of eigenangles in previous segment
    xmtrfields = Vector{T}(undef, 0)  # fields generated by transmitter
    previous_xmtrfields = similar(xmtrfields)  # fields saved from previous segment
    rcvrfields = similar(xmtrfields)  # fields at receiver location

    Xidx = 1
    for segmentidx = 1:numsegments
        wvg = waveguide[segmentidx]
        wavefields = wavefields_vec[segmentidx]
        eas = eigenangles(wavefields)
        current_eacount = length(eas)

        # Identify distance at beginning of segment
        segment_start = wvg.distance
        if segmentidx == 1
            @assert segment_start == 0 "First waveguide segment should have distance 0"
        end
        maxX < segment_start && break  # no farther X; break

        # Identify distance at end of segment
        if segmentidx < numsegments
            segment_end = waveguide[segmentidx+1].distance
        else
            segment_end = typemax(typeof(segment_start))
        end

        # xmtrfields is for `Hy`
        resize!(xmtrfields, current_eacount)
        resize!(rcvrfields, current_eacount)
        if segmentidx > 1
            adjwavefields = adjwavefields_vec[segmentidx]
            prevwavefields = wavefields_vec[segmentidx-1]
            conversioncoeffs = modeconversion(prevwavefields, wavefields, adjwavefields)
        end

        # Calculate the mode terms (height gains and excitation factors) up to
        # the current segment
        for n = 1:current_eacount
            # `xmtrterm` includes excitation factor and transmitter height gain
            # `rcvrterm` is receiver height gain only
            modeequation = PhysicalModeEquation(eas[n], frequency, wvg)
            xmtrterm, rcvrterm, efc = modeterms(modeequation, tx, rx)
            if segmentidx == 1
                # Transmitter exists only in the transmitter slab (obviously)
                xmtrfields[n] = xmtrterm
            else
                # Otherwise, mode conversion of transmitted fields
                xmtrfields_sum = zero(eltype(xmtrfields))
                for m = 1:previous_eacount
                    xmtrfields_sum += previous_xmtrfields[m]*conversioncoeffs[n,m]  # TODO: swap conversioncoeffs n,m?
                end
                xmtrfields[n] = xmtrfields_sum
            end

            rcvrfields[n] = xmtrfields[n]*rcvrterm
        end

        # Calculate E at each distance in the current waveguide segment
        while X[Xidx] < segment_end
            x = X[Xidx] - segment_start
            factor = Q/sqrt(abs(sin(X[Xidx]/EARTHRADIUS)))

            totalfield = zero(T)
            for n = 1:current_eacount
                S₀ = referencetoground(eas[n].sinθ)
                totalfield += rcvrfields[n]*cis(-k*x*(S₀ - 1))*factor
            end

            E[Xidx] = totalfield
            Xidx += 1
            Xidx > Xlength && break
        end

        # If we've reached the end of the current segment and there are more
        # segments, prepare for next segment
        if segmentidx < numsegments
            # End of current slab
            x = segment_end - segment_start

            resize!(previous_xmtrfields, current_eacount)
            for n = 1:current_eacount
                S₀ = referencetoground(eas[n].sinθ)

                # Excitation factors at end of slab. LWPC uses `Hy`
                xmtrfields[n] *= cis(-k*x*(S₀ - 1))
                previous_xmtrfields[n] = xmtrfields[n]
            end
            previous_eacount = current_eacount
        end
    end

    return E
end


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


########
# Utility functions

function unwrap!(x)
	v = first(x)
	@inbounds for k in eachindex(x)
		x[k] = v = v + rem2pi(x[k]-v, RoundNearest)
	end
	return x
end

# see, e.g. PS71 pg 11
function referencetoground(ea::EigenAngle)
    return EigenAngle(asin(ea.sinθ/sqrt(1 - 2/EARTHRADIUS*CURVATURE_HEIGHT)))
end
referencetoground(x::Number) = x/sqrt(1 - 2/EARTHRADIUS*CURVATURE_HEIGHT)

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
