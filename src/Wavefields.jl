#==
Functions related to calculating electromagnetic field components at any height
within the ionosphere.
==#

"""
    Wavefields{H<:AbstractRange}

Struct containing wavefield components for each mode of `eas` at `heights`.

The six electromagnetic field components are stored in the `v` field of the struct and
accessed as a `Matrix` of `SVector{6,ComplexF64}` corresponding to `[height,mode]`.
"""
struct Wavefields{H<:AbstractRange}
    v::Matrix{SVector{6,ComplexF64}}  # v[height,mode]  ; type is derived from typeof T, which is from M
    heights::H
    eas::Vector{EigenAngle}
end

"""
    Wavefields(eas::Vector{EigenAngle}, heights)

A constructor for initializing `Wavefields` with an appropriately sized `undef` Matrix
given eigenangles `eas` and `heights`.
"""
function Wavefields(heights, eas::Vector{EigenAngle})
    v = Matrix{SVector{6,ComplexF64}}(undef, length(heights), length(eas))
    Wavefields(v, heights, eas)
end

Base.size(A::Wavefields) = size(A.v)
Base.size(A::Wavefields, I...) = size(A.v, I...)
Base.getindex(A::Wavefields, i) = A.v[i]
Base.getindex(A::Wavefields, I...) = A.v[I...]
Base.setindex!(A::Wavefields, x, i) = (A.v[i] = x)
Base.similar(A::Wavefields) = Wavefields(A.heights, A.eas)
Base.copy(A::Wavefields) = Wavefields(copy(A.v), copy(A.heights), copy(A.eas))
Base.view(A::Wavefields, I...) = view(A.v, I...)

function (==)(A::Wavefields, B::Wavefields)
    A.heights == B.heights || return false
    A.eas == B.eas || return false
    A.v == B.v || return false
    return true
end

heights(A::Wavefields) = A.heights
numheights(A::Wavefields) = length(A.heights)
eigenangles(A::Wavefields) = A.eas
nummodes(A::Wavefields) = length(A.eas)

function Base.isvalid(A::Wavefields)
    vzs, veas = size(A.v)
    zlen = numheights(A)
    vzs == zlen || return false
    ealen = nummodes(A)
    veas == ealen || return false
    return true
end


"""
    WavefieldIntegrationParams{F,G,H}(z, bottomz, ortho_scalar, e1_scalar, e2_scalar, ea,
        frequency, bfield, species, lwmsparams)

Parameters passed to Pitteway integration of wavefields.

    - `z::Float64`
    - `bottomz::Float64`
    - `ortho_scalar::Complex{Float64}`
    - `e1_scalar::Float64`
    - `e2_scalar::Float64`
    - `ea::EigenAngle`
    - `frequency::Frequency`
    - `bfield::BField`
    - `species::Species{F,G}`
    - `lwmsparams::LWMSParams{H}`
"""
struct WavefieldIntegrationParams{F,G,H}
    z::Float64
    bottomz::Float64
    ortho_scalar::ComplexF64
    e1_scalar::Float64
    e2_scalar::Float64
    ea::EigenAngle
    frequency::Frequency
    bfield::BField
    species::Species{F,G}
    lwmsparams::LWMSParams{H}
end

"""
    WavefieldIntegrationParams(topheight, ea, frequency, bfield, species, lwmsparams)

Initialize a `WavefieldIntegrationParams` for downward Pitteway scaled integration.

Automatically set values are:

    - `bottomz = BOTTOMHEIGHT`
    - `ortho_scalar = zero(ComplexF64)`
    - `e1_scalar = one(Float64)`
    - `e2_scalar = one(Float64)`
"""
function WavefieldIntegrationParams(topheight, ea, frequency, bfield, species::Species{F,G},
    lwmsparams::LWMSParams{H}) where {F,G,H}
    return WavefieldIntegrationParams{F,G,H}(topheight, BOTTOMHEIGHT, zero(ComplexF64),
        one(Float64), one(Float64), ea, frequency, bfield, species, lwmsparams)
end

"""
    ScaleRecord(z, e, ortho_scalar, e1_scalar, e2_scalar)

Struct used for saving wavefield scaling information during Pitteway integration of
wavefields.
"""
struct ScaleRecord
    z::Float64
    e::SMatrix{4,2,Complex{Float64},8}
    ortho_scalar::Complex{Float64}
    e1_scalar::Float64
    e2_scalar::Float64
end

"""
    bookerwavefields(T::TMatrix)

Calculate the initial wavefields vector ``[Ex₁ Ex₂
                                           Ey₁ Ey₂
                                           Hx₁ Hx₂
                                           Hy₁ Hy₂]``
for the two upgoing wavefields where subscript `1` is the evanescent wave and
`2` is the travelling wave.

The wavefields are always of type `Complex{Float64}`.

This function solves the equation ``Te = qe`` equivalent to the eigenvector
problem ``(T- qI)e = 0``. First, the Booker Quartic is solved for the roots `q`,
and they are sorted so that the roots associated with the two upgoing waves are
selected, where eigenvalue ``q₁`` corresponds to the evanescent wave and ``q₂``
the travelling wave. Then `e` is solved as the eigenvectors for the two `q`s. An
analytical solution is used where `e[2,:] = 1`.
"""
function bookerwavefields(T::TMatrix)
    q, B = bookerquartic(T)
    sortquarticroots!(q)
    return bookerwavefields(T, q)
end

function bookerwavefields(ea::EigenAngle, M)
    q, B = bookerquartic(ea, M)
    sortquarticroots!(q)
    T = tmatrix(ea, M)
    return bookerwavefields(T, q)
end

function bookerwavefields(T::TMatrix, q)
    # Precompute
    T14T41 = T[1,4]*T[4,1]
    T14T42 = T[1,4]*T[4,2]
    T12T41 = T[1,2]*T[4,1]

    e = MArray{Tuple{4,2},ComplexF64,2,8}(undef)
    @inbounds for i = 1:2
        d = T14T41 - (T[1,1] - q[i])*(T[4,4] - q[i])

        e[1,i] = (T[1,2]*(T[4,4] - q[i]) - T14T42)/d
        e[2,i] = 1
        e[3,i] = q[i]
        e[4,i] = (-T12T41 + T[4,2]*(T[1,1] - q[i]))/d
    end

    # By returning as SArray instead of MArray, the MArray doesn't get hit by GC
    return SArray(e)
end

function bookerwavefields(T::TMatrix, dT, ::Dθ)
    q, B = bookerquartic(T)
    sortquarticroots!(q)
    dq = bookerquartic(T, dT, q, B)
    return bookerwavefields(T, dT, q, dq)
end

function bookerwavefields(ea::EigenAngle, M, ::Dθ)
    q, B = bookerquartic(ea, M)
    sortquarticroots!(q)
    dq = bookerquartic(ea, M, q, B)
    T = tmatrix(ea, M)
    dT = tmatrix(ea, M, Dθ())

    return bookerwavefields(T, dT, q, dq)
end

function bookerwavefields(T::TMatrix, dT, q, dq)
    # Precompute
    T14T41 = T[1,4]*T[4,1]
    T14T42 = T[1,4]*T[4,2]
    T12T41 = T[1,2]*T[4,1]

    e = MArray{Tuple{4,2},ComplexF64,2,8}(undef)
    de = similar(e)
    @inbounds for i = 1:2
        d = T14T41 - (T[1,1] - q[i])*(T[4,4] - q[i])
        d2 = abs2(d)
        dd = dT[1,4]*T[4,1] - T[1,1]*dT[4,4] - dT[1,1]*T[4,4] +
            q[i]*(dT[1,1] + dT[4,4] - 2*dq[i]) + dq[i]*(T[1,1] + T[4,4])

        e[1,i] = (T[1,2]*(T[4,4] - q[i]) - T14T42)/d
        e[2,i] = 1
        e[3,i] = q[i]
        e[4,i] = (-T12T41 + T[4,2]*(T[1,1] - q[i]))/d

        de1num = (T[1,2]*(dT[4,4] - dq[i]) + dT[1,2]*(T[4,4] - q[i]) - dT[1,4]*T[4,2])
        de[1,i] = (de1num*d - e[1,i]*d*dd)/d2
        de[2,i] = 0
        de[3,i] = dq[i]
        de4num = -dT[1,2]*T[4,1] + T[4,2]*(dT[1,1] - dq[i])
        de[4,i] = (de4num*d - e[4,i]*d*dd)/d2
    end

    # By returning as SArray instead of MArray, the MArray doesn't get hit by GC
    return SArray(e), SArray(de)
end

"""
    dedz(e, k, T::Tmatrix)

Calculates derivative of field components vector ``de/dz = -i k T e``.
"""
dedz(e, k, T::TMatrix) = -1im*k*(T*e)  # `(T*e)` uses specialized TMatrix math

"""
    dedz(e, p, z)

Calculates derivative of field components vector `e` at height `z`.

The parameters `p` should be a `WavefieldIntegrationParams`.

This function internally calls [`susceptibility`](@ref) and [`tmatrix`](@ref) and is typically
used by [`integratewavefields`](@ref).
"""
function dedz(e, p, z)
    @unpack ea, frequency, bfield, species, lwmsparams = p

    M = susceptibility(z, frequency, bfield, species, params=lwmsparams)
    T = tmatrix(ea, M)

    return dedz(e, frequency.k, T)
end

"""
    scalingcondition(e, z, int)

Return `true` if wavefields should be scaled, otherwise `false` for wavefields `e` at height
`z` and integrator `int`.

Specifically, if any component of `real(e)` or `imag(e)` are `>= 1`, return
`true`. In addition, force scale at `z == bottomz` to ensure initial upgoing wave
is unit amplitude.
"""
scalingcondition(e, z, int) = any(x->(real(x) >= 1 || imag(x) >= 1), e) || z == int.p.bottomz

"""
    scale!(integrator)

Apply wavefield scaling with [`scalewavefields`](@ref) to the integrator.
"""
function scale!(integrator)
    new_e, new_orthos, new_e1s, new_e2s = scalewavefields(integrator.u)

    # Last set of scaling values
    @unpack bottomz, ea, frequency, bfield, species, lwmsparams = integrator.p

    #==
    NOTE: `integrator.t` is the "time" of the _proposed_ step. Therefore,
    integrator.t` might equal `0.0`, for example, before it's actually gotten
    to the bottom. `integrator.prevt` is the last `t` on the "left"
    side of the `integrator`, which covers the local interval [`tprev`, `t`].
    The "condition" is met at `integrator.t` and `integrator.t` is the time at
    which the affect occurs.
    However, it is not guaranteed that each (`tprev`, `t`) directly abuts the
    next `tprev`, `t`).
    ==#

    integrator.p = WavefieldIntegrationParams(integrator.t,
                                              bottomz,
                                              new_orthos,
                                              new_e1s, new_e2s,
                                              ea, frequency, bfield, species, lwmsparams)

    integrator.u = new_e

    return nothing
end

"""
    scalewavefields(e1, e2)

Returns orthonormalized vectors `e1` and `e2`, as well as the scaling terms `a`,
`e1_scale_val`, and `e2_scale_val` applied to the original vectors.

First applies Gram-Schmidt orthogonalization and then scales the vectors so they
each have length 1, i.e. `norm(e1) == norm(e2) == 1`. This is the technique
suggested by [^Pitteway1965] to counter numerical swamping during integration of
wavefields.

# References

[^Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields,
reflexion coefficients and polarizations for long radio waves in the lower
ionosphere. I.,” Phil. Trans. R. Soc. Lond. A, vol. 257, no. 1079,
pp. 219–241, Mar. 1965, doi: 10.1098/rsta.1965.0004.
"""
function scalewavefields(e1, e2)
    # Orthogonalize vectors `e1` and `e2` (Gram-Schmidt process)
    # `dot` for complex vectors automatically conjugates first vector
    e1_dot_e1 = real(dot(e1, e1))  # == sum(abs2.(e1)), `imag(dot(e1,e1)) == 0`
    a = dot(e1, e2)/e1_dot_e1  # purposefully unsigned
    e2 -= a*e1

    # Normalize `e1` and `e2`
    e1_scale_val = inv(sqrt(e1_dot_e1))
    e2_scale_val = inv(norm(e2))  # == 1/sqrt(dot(e2,e2))
    e1 *= e1_scale_val
    e2 *= e2_scale_val  # == normalize(e2)

    return e1, e2, a, e1_scale_val, e2_scale_val
end

"""
    scalewavefields(e)

!!! note

    This function only applies scaling to the first 2 columns of `e`.
"""
function scalewavefields(e)
    e1, e2, a, e1_scale_val, e2_scale_val = scalewavefields(e[:,1], e[:,2])

    if size(e, 2) == 2
        e = hcat(e1, e2)
    else
        e = hcat(e1, e2, e[:,3:end])
    end

    return e, a, e1_scale_val, e2_scale_val
end

"""
    unscalewavefields!(e, saved_values::SavedValues)

Unscale the integrated wavefields `e` in place.

Assumes fields have been scaled by [`scalewavefields`](@ref) during integration.

See also [`unscalewavefields`](@ref) for additional details.
"""
function unscalewavefields!(e, saved_values::SavedValues)
    zs = saved_values.t
    records = saved_values.saveval

    # Initialize the "reference" scaling altitude
    ref_z = last(records).z

    # Usually `osum = 0`, `prod_e1 = 1`, and `prod_e2 = 1` at initialization,
    # but we set up the fields at the ground (`last(records)`) outside the loop
    osum = last(records).ortho_scalar
    prod_e1 = last(records).e1_scalar
    prod_e2 = last(records).e2_scalar

    # Unscaling we go from the bottom up
    @inbounds for i in reverse(eachindex(e))

        # Unpack variables
        record_z = records[i].z
        scaled_e = records[i].e
        ortho_scalar = records[i].ortho_scalar
        e1_scalar = records[i].e1_scalar
        e2_scalar = records[i].e2_scalar

        z = zs[i]

        # Only update ref_z when there is both:
        # 1) we have reached the height where a new ref_z should go into effect
        # 2) there is a new ref_z
        if (z > record_z) & (record_z > ref_z)
            ref_z = record_z

            osum *= e1_scalar/e2_scalar
            osum += ortho_scalar
            prod_e1 *= e1_scalar
            prod_e2 *= e2_scalar
        end

        if i == lastindex(e)
            # Bottom doesn't require correction
            e[i] = scaled_e
        elseif z > ref_z
            # From the bottom, the first correction may not need to be applied
            # until some higher altitude
            e2 = (scaled_e[:,2] - osum*scaled_e[:,1])*prod_e2
            e1 = scaled_e[:,1]*prod_e1
            e[i] = hcat(e1,e2)
        else
            e[i] = scaled_e
        end
    end

    return nothing
end

"""
    unscalewavefields(saved_values::SavedValues)

Return the unscaled integrated wavefields originally scaled by [`scalewavefields`](@ref).

The fields are not "unscaled" so much as further scaled so that at all heights the fields
are scaled similarly to the cumulative scaling that applied to the fields at the bottom.

The bottom level does not get unscaled. We reference the higher levels to the bottom. The
level above the bottom level needs to be additionally scaled by the amount that was applied
to originally get from this level down to the bottom level. The next level up (2 above the
bottom level) needs to be scaled by the amount applied to the next level and then the bottom
level, i.e. we keep track of a cumulative correction on the way back up.
"""
function unscalewavefields(saved_values::SavedValues)
    e = Vector{SArray{Tuple{4,2},Complex{Float64},2,8}}(undef, length(saved_values.saveval))

    unscalewavefields!(e, saved_values)

    return e
end

"""
    save_values(u, t, integrator)

Return a `ScaleRecord` from `u`, `t`, and `integrator`.

Used by SavingCallback in [`integratewavefields`](@ref).
"""
save_values(u, t, integrator) = ScaleRecord(integrator.p.z,
                                            u,
                                            integrator.p.ortho_scalar,
                                            integrator.p.e1_scalar,
                                            integrator.p.e2_scalar)

"""
    integratewavefields(zs, ea::EigenAngle, frequency::Frequency, bfield::BField,
        species::Species; params=LWMSParams())

Return wavefields `e` at `zs` by downward integration over `zs`.
"""
function integratewavefields(zs, ea::EigenAngle, frequency::Frequency, bfield::BField,
    species::Species; params=LWMSParams())
    # TODO: version that updates output `e` in place

    issorted(zs, rev=true) || throw(ArgumentError("zs should go from top to bottom"))

    # Initial conditions
    Mtop = susceptibility(first(zs), frequency, bfield, species, params=params)
    Ttop = tmatrix(ea, Mtop)
    e0 = bookerwavefields(Ttop)

    # Works best if `e0` fields are normalized at top height
    # Don't want to orthogonalize `e0[:,2]` here because it won't be undone!
    e0 = hcat(normalize(e0[:,1]), normalize(e0[:,2]))

    # saved_positions=(true, true) because we discontinuously modify `u`.
    # This is independent of `saveat` and `save_everystep`
    cb = DiscreteCallback(scalingcondition, scale!, save_positions=(true, true))

    saved_values = SavedValues(Float64, ScaleRecord)

    scb = SavingCallback(save_values, saved_values, save_everystep=false, saveat=zs, tdir=-1)

    p = WavefieldIntegrationParams(params.topheight, ea, frequency, bfield, species, params)

    # WARNING: Without `lazy=false` (necessary since we're using DiscreteCallback) don't
    # use continuous solution output! Also, we're only saving at zs.
    prob = ODEProblem{false}(dedz, e0, (first(zs), last(zs)), p)
    sol = solve(prob, Vern9(lazy=false), callback=CallbackSet(cb, scb),
                save_everystep=false, save_start=false, save_end=false,
                atol=1e-8, rtol=1e-8)

    e = unscalewavefields(saved_values)

    return e
end

"""
    vacuumreflectioncoeffs(ea, e1, e2)

Return ionosphere reflection coefficient matrix from upgoing wave fields `e`.

Integrating for one set of horizontal field components ``e = (Ex, -Ey, Z₀Hx, Z₀Hy)ᵀ``
can be separated into an upgoing and downgoing wave, each of which is generally
elliptically polarized. One might assume that the ratio of the amplitudes of
these two waves would give a reflection coefficient, and it does, except the
coefficient would only apply for an incident wave of that particular elliptical
polarization. However, the first set of fields can be linearly combined with
a second independent solution for the fields, which will generally have a
different elliptical polarization than the first. Two linear combinations of the
two sets of fields are formed with unit amplitude, linearly polarized
incident waves. The reflected waves then give the components ``R₁₁``, ``R₂₁`` or
``R₁₂``, ``R₂₂`` for the incident wave in the plane of incidence and
perpendicular to it, respectively [^Budden1988] (pg 552).

The process for determining the reflection coefficient requires resolving the
two sets of fields `e1` and `e2` into the four linearly polarized vacuum
modes. The layer of vacuum can be assumed to be so thin that it does not affect
the fields. There will be two upgoing waves, one of which has ``E``, and the
other ``H`` in the plane of incidence, and two downgoing waves, with ``E`` and
``H`` in the plane of incidence. If ``f₁, f₂, f₃, f₄`` are the complex
amplitudes of the four component waves, then in matrix notation ``e = Sᵥ f``.

For `e1` and `e2`, we can find the corresponding vectors `f1` and `f2` by
``f1 = Sᵥ⁻¹ e1``, ``f2 = Sᵥ⁻¹ e2`` where the two column vectors are partitioned
such that ``f1 = (u1, d1)ᵀ`` and ``f2 = (u2, d2)ᵀ`` for upgoing and downgoing
2-element vectors `u` and `d`. From the definition of the reflection coefficient
`R`, ``d = Ru``. Letting ``U = (u1, u2)``, ``D = (d1, d2)``, then ``D = RU`` and
the reflection coefficient is ``R = DU¹``. Because the reflection coefficient
matrix is a ratio of fields, either `e1` and/or `e2` can be independently
multiplied by an arbitrary constant and the value of `R` is unaffected.

For additional details, see [^Budden1988], chapter 18, section 7.

# References

[^Budden1988] K. G. Budden, The propagation of radio waves: the theory of radio
waves of low power in the ionosphere and magnetosphere, First paperback edition.
New York: Cambridge University Press, 1988.
"""
vacuumreflectioncoeffs(ea::EigenAngle, e1, e2) = bookerreflection(ea, hcat(e1, e2))

"""
    vacuumreflectioncoeffs(ea, e)

Calculate `vacuumreflectioncoeffs` where `e` can be decomposed into `e[:,1]` and
`e[:,2]`.
"""
vacuumreflectioncoeffs(ea, e) = bookerreflection(ea, e)

"""
    wavefieldboundary(R, Rg, e1, e2)

Calculate coefficients `b1`, `b2` required to sum `e1`, `e2` for total wavefield
as ``e = b1*e1 + b2*e2``.

This process is used by LWPC ("wf_bndy.FOR") and is similar to the calculation of excitation
factor at the ground because it makes use of the mode condition.
"""
function boundaryscalars(R, Rg, e1, e2, isotropic=false)
    # This is similar to excitation factor calculation, using the waveguide mode condition

    ex1, ey1, hx1, hy1 = e1[1], -e1[2], e1[3], e1[4]
    ex2, ey2, hx2, hy2 = e2[1], -e2[2], e2[3], e2[4]

    # TODO modify for flat earth case (imag < -10)
    # NOTE: these assume `R`s and `e`s at the ground, see e.g. "wf_rbars.for"
    hy0 = 1  # == (1 + Rg[1,1])/(1 + Rg[1,1])
    ey0 = 1  # == (1 + Rg[2,2])/(1 + Rg[2,2])

    abparal = 1 - Rg[1,1]*R[1,1]
    abperp = 1 - Rg[2,2]*R[2,2]

    if isotropic
        if abs2(abperp) < abs2(abparal)
            b1 = zero(ey2)
            b2 = ey0/ey2
            # hyg = 0
        else
            b1 = hy0/hy1
            b2 = zero(hy1)
            # eyg = 0
        end
    else
        # polarization ratio Ey/Hy (often `f` or `fofr` in papers)
        if abs2(abparal) < abs2(abperp)
            pol = ((1 + Rg[2,2])*Rg[1,1]*R[2,1])/((1 + Rg[1,1])*abperp)
        else
            pol = ((1 + Rg[2,2])*abparal)/((1 + Rg[1,1])*Rg[2,2]*R[1,2])
        end
        a = (-ey1 + pol*hy1)/(ey2 - pol*hy2)
        hysum = hy1 + a*hy2
        b1 = hy0/hysum
        b2 = b1*a
    end

    return b1, b2
end

boundaryscalars(R, Rg, e, isotropic::Bool) = boundaryscalars(R, Rg, e[:,1], e[:,2], isotropic)

"""
    fieldstrengths()

Calculate `Ex`, `Ey`, `Ez`, `Hx`, `Hy`, `Hz` wavefields by full wave integration at heights
`zs`.

Scales wavefields for waveguide boundary conditions.
"""
function fieldstrengths!(EH, zs, ea::EigenAngle, frequency::Frequency, bfield::BField,
    species::Species, ground::Ground; params=LWMSParams())

    @assert length(EH) == length(zs)
    zs[end] == 0 || @warn "Waveguide math assumes fields and reflection coefficients are
        calculated at the ground (`z = 0`)."

    e = integratewavefields(zs, ea, frequency, bfield, species, params=params)
    R = vacuumreflectioncoeffs(ea, e[end])
    Rg = fresnelreflection(ea, ground, frequency)
    b1, b2 = boundaryscalars(R, Rg, e[end], isisotropic(bfield))

    S = ea.sinθ
    @inbounds for i in eachindex(EH)
        # Scale to waveguide boundary conditions
        w = e[i][:,1]*b1 + e[i][:,2]*b2

        M = susceptibility(zs[i], frequency, bfield, species, params=params)
        ez = -(w[4]*S + M[3,1]*w[1] - M[3,2]*w[2])/(1 + M[3,3])
        hz = -w[2]*S
        EH[i] = SVector(w[1], -w[2], ez, w[3], w[4], hz)
    end

    return nothing
end

function calculate_wavefields!(wavefields, adjoint_wavefields, frequency, waveguide,
    adjoint_waveguide; params=LWMSParams())

    @unpack bfield, species, ground = waveguide

    ground == adjoint_waveguide.ground || throw(ArgumentError("waveguide and adjoint_waveguide should have same Ground"))
    species == adjoint_waveguide.species || throw(ArgumentError("waveguide and adjoint_waveguide should have same Species"))

    adjoint_bfield = adjoint_waveguide.bfield

    zs = heights(wavefields)
    modes = eigenangles(wavefields)
    adjoint_modes = eigenangles(adjoint_wavefields)
    modes == adjoint_modes || @warn "Full mode conversion physics assumes adjoint
        wavefields are calculated for eigenangles of the original waveguide."

    for m in eachindex(modes)
        fieldstrengths!(view(wavefields,:,m), zs, modes[m], frequency, bfield, species,
                        ground, params=params)
        fieldstrengths!(view(adjoint_wavefields,:,m), zs, adjoint_modes[m], frequency,
                        adjoint_bfield, species, ground, params=params)
    end

    return nothing
end
