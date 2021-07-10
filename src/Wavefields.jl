#==
Functions related to calculating electromagnetic field components at any height within the
ionosphere.

This includes calculating wavefields vector from the Booker quartic and Pitteway
integration of wavefields.
==#

"""
    Wavefields{H<:AbstractRange}

Struct containing wavefield components for each mode of `eas` at `heights`.

The six electromagnetic field components are stored in the `v` field of the struct and
accessed as a `Matrix` of `SVector{6,ComplexF64}` corresponding to `[height, mode]`.
"""
struct Wavefields{H<:AbstractRange}
    v::Matrix{SVector{6,ComplexF64}}  # v[height,mode]  ; type is derived from typeof T, which is from M
    heights::H
    eas::Vector{EigenAngle}
end

"""
    Wavefields(heights, eas::Vector{EigenAngle})

A constructor for initializing `Wavefields` with an appropriately sized `undef` Matrix
given eigenangles `eas` and `heights`.
"""
function Wavefields(heights, eas::Vector{EigenAngle})
    v = Matrix{SVector{6,ComplexF64}}(undef, length(heights), length(eas))
    return Wavefields(v, heights, eas)
end

Base.size(A::Wavefields) = size(A.v)
Base.size(A::Wavefields, I...) = size(A.v, I...)
Base.getindex(A::Wavefields, i) = A.v[i]
Base.getindex(A::Wavefields, I...) = A.v[I...]
Base.setindex!(A::Wavefields, x, i) = (A.v[i] = x)
Base.similar(A::Wavefields) = Wavefields(A.heights, A.eas)
Base.copy(A::Wavefields) = Wavefields(copy(A.v), copy(A.heights), copy(A.eas))
Base.view(A::Wavefields, I...) = view(A.v, I...)

"""
    ==(A::Wavefields, B::Wavefields)

Check for equality `==` between `A.v` and `B.v`, `heights`, and `eas`.
"""
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
    WavefieldIntegrationParams{S,T,T2,H}

Parameters passed to Pitteway integration of wavefields [Pitteway1965].

# Fields

- `z::Float64`: height in meters.
- `bottomz::Float64`: bottom height of integration in meters.
- `ortho_scalar::Complex{Float64}`: scalar for Gram-Schmidt orthogonalization.
- `e1_scalar::Float64`: scalar for wavefield vector 1.
- `e2_scalar::Float64`: scalar for wavefield vector 2.
- `ea::EigenAngle`: wavefield eigenangle.
- `frequency::Frequency`: electromagentic wave frequency in Hz.
- `bfield::BField`: Earth's magnetic field in Tesla.
- `species::S`: species in the ionosphere.
- `params::LMPParams{T,T2,H}`: module-wide parameters.
"""
struct WavefieldIntegrationParams{S,T,T2,H}
    z::Float64
    bottomz::Float64
    ortho_scalar::ComplexF64
    e1_scalar::Float64
    e2_scalar::Float64
    ea::EigenAngle
    frequency::Frequency
    bfield::BField
    species::S
    params::LMPParams{T,T2,H}
end

"""
    WavefieldIntegrationParams(topheight, ea, frequency, bfield, species, params)

Initialize a `WavefieldIntegrationParams` for downward Pitteway scaled integration
[Pitteway1965].

Automatically set values are:

- `bottomz = BOTTOMHEIGHT`
- `ortho_scalar = zero(ComplexF64)`
- `e1_scalar = one(Float64)`
- `e2_scalar = one(Float64)`
"""
function WavefieldIntegrationParams(topheight, ea, frequency, bfield, species::S,
    params::LMPParams{T,T2,H}) where {S,T,T2,H}
    return WavefieldIntegrationParams{S,T,T2,H}(topheight, BOTTOMHEIGHT, zero(ComplexF64),
        1.0, 1.0, ea, frequency, bfield, species, params)
end

"""
    ScaleRecord

Struct for saving wavefield scaling information used in a callback during Pitteway
integration of wavefields [Pitteway1965].

# Fields

- `z::Float64`: current height in meters.
- `e::SMatrix{4,2,Complex{Float64},8}`: wavefield matrix.
- `ortho_scalar::Complex{Float64}`: scalar for Gram-Schmidt orthogonalization.
- `e1_scalar::Float64`: scalar for wavefield vector 1.
- `e2_scalar::Float64`: scalar for wavefield vector 2.
"""
struct ScaleRecord
    z::Float64
    e::SMatrix{4,2,Complex{Float64},8}
    ortho_scalar::Complex{Float64}
    e1_scalar::Float64
    e2_scalar::Float64
end

"""
    dedz(e, k, T::Tmatrix)

Compute ``de/dz = -ikTe`` where ``e = (Ex, -Ey, Z₀Hx, Z₀Hy)ᵀ``.
"""
dedz(e, k, T::TMatrix) = -1im*k*(T*e)  # `(T*e)` uses specialized TMatrix math

"""
    dedz(e, p, z)

Compute derivative of field components vector `e` at height `z`.

The parameter `p` should be a `WavefieldIntegrationParams`.
"""
function dedz(e, p, z)
    @unpack ea, frequency, bfield, species, params = p

    M = susceptibility(z, frequency, bfield, species, params=params)
    T = tmatrix(ea, M)

    return dedz(e, frequency.k, T)
end

"""
    scalingcondition(e, z, int)

Return `true` if wavefields should be scaled, otherwise `false` for wavefields `e` at height
`z` and integrator `int`.

Wavefields should be scaled if any component of `real(e)` or `imag(e)` are `>= 1`. In
addition, force scaling at `z <= bottomz` to ensure initial upgoing wave is unit amplitude.
"""
scalingcondition(e, z, int) = any(x->(real(x) >= 1 || imag(x) >= 1), e) || z <= int.p.bottomz
# scalingcondition(e, z, int) = any(x -> abs2(x) >= 1, e) || z <= int.p.bottomz

"""
    scale!(integrator)

Apply wavefield scaling with [`scalewavefields`](@ref) to the integrator.
"""
function scale!(integrator)
    new_e, new_orthos, new_e1s, new_e2s = scalewavefields(integrator.u)

    # Last set of scaling values
    @unpack bottomz, ea, frequency, bfield, species, params = integrator.p

    #==
    NOTE: `integrator.t` is the "time" of the _proposed_ step. Therefore, integrator.t`
    might equal `0.0`, for example, before it's actually gotten to the bottom.
    `integrator.prevt` is the last `t` on the "left" side of the `integrator`, which covers
    the local interval [`tprev`, `t`]. The "condition" is met at `integrator.t` and
    `integrator.t` is the time at which the affect occurs.
    However, it is not guaranteed that each (`tprev`, `t`) directly abuts the next (`tprev`,
    `t`).
    ==#

    integrator.p = WavefieldIntegrationParams(integrator.t,
                                              bottomz,
                                              new_orthos,
                                              new_e1s, new_e2s,
                                              ea, frequency, bfield, species, params)

    integrator.u = new_e

    return nothing
end

"""
    scalewavefields(e1, e2)
    scalewavefields(e::SMatrix{4,2})

Orthonormalize the vectors `e1` and `e2` or the columns of `e`, and also return the scaling
terms `a`, `e1_scale_val`, and `e2_scale_val` applied to the original vectors.

This first applies Gram-Schmidt orthogonalization and then scales the vectors so they each
have length 1, i.e. `norm(e1) == norm(e2) == 1`. This is the technique suggested by
[Pitteway1965] to counter numerical swamping during integration of wavefields.

# References

[Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields, reflexion
    coefficients and polarizations for long radio waves in the lower ionosphere. I.,” Phil.
    Trans. R. Soc. Lond. A, vol. 257, no. 1079, pp. 219–241, Mar. 1965,
    doi: 10.1098/rsta.1965.0004.
"""
scalewavefields

function scalewavefields(e1, e2)
    # Orthogonalize vectors `e1` and `e2` (Gram-Schmidt process)
    # `dot` for complex vectors automatically conjugates first vector
    e1_dot_e1 = real(dot(e1, e1))  # == sum(abs2.(e1)), imag(dot(e1,e1)) == 0
    a = dot(e1, e2)/e1_dot_e1  # purposefully unsigned
    e2 -= a*e1

    # Normalize `e1` and `e2`
    e1_scale_val = 1/(sqrt(e1_dot_e1))
    e2_scale_val = 1/(norm(e2))  # == 1/sqrt(dot(e2,e2))
    e1 *= e1_scale_val
    e2 *= e2_scale_val  # == normalize(e2)

    return e1, e2, a, e1_scale_val, e2_scale_val
end

function scalewavefields(e::SMatrix{4,2})
    e1, e2, a, e1_scale_val, e2_scale_val = scalewavefields(e[:,1], e[:,2])

    e = hcat(e1, e2)

    return e, a, e1_scale_val, e2_scale_val
end

"""
    unscalewavefields!(e, saved_values::SavedValues)

Unscale the integrated wavefields, a vector of fields at each integration step `e`, in place.

See also: [`unscalewavefields`](@ref)
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

The bottom level does not get unscaled. We reference the higher levels to the bottom. The
level above the bottom level is additionally scaled by the amount that was applied to
originally get from this level down to the bottom level. The next level up (2 above the
bottom level) is scaled by the amount applied to the next level and then the bottom
level, i.e. we keep track of a cumulative correction on the way back up.

Assumes fields have been scaled by [`scalewavefields`](@ref) during integration.

See also: [`unscalewavefields!`](@ref)
"""
function unscalewavefields(saved_values::SavedValues)
    e = Vector{SArray{Tuple{4,2},Complex{Float64},2,8}}(undef, length(saved_values.saveval))

    unscalewavefields!(e, saved_values)

    return e
end

"""
    savevalues(u, t, integrator)

Return a `ScaleRecord` from `u`, `t`, and `integrator`.
"""
savevalues(u, t, integrator) = ScaleRecord(integrator.p.z,
                                           u,
                                           integrator.p.ortho_scalar,
                                           integrator.p.e1_scalar,
                                           integrator.p.e2_scalar)

"""
    integratewavefields(zs, ea::EigenAngle, frequency::Frequency, bfield::BField,
        species; params=LMPParams())

Compute wavefields vectors `e` at `zs` by downward integration over heights `zs`.

`params.wavefieldintegrationparams` is used by this function rather than
`params.integrationparams`.
"""
function integratewavefields(zs, ea::EigenAngle, frequency::Frequency, bfield::BField,
    species; params=LMPParams(), unscale=true)
    # TODO: version that updates output `e` in place

    issorted(zs; rev=true) ||
        throw(ArgumentError("`zs` should go from top to bottom of the ionosphere."))

    @unpack wavefieldintegrationparams = params
    @unpack tolerance, solver, dt, force_dtmin, maxiters = wavefieldintegrationparams

    # Initial conditions
    Mtop = susceptibility(first(zs), frequency, bfield, species; params=params)
    Ttop = tmatrix(ea, Mtop)
    e0 = bookerwavefields(Ttop)

    # Works best if `e0` fields are normalized at top height
    # Don't want to orthogonalize `e0[:,2]` here because it won't be undone!
    e0 = hcat(normalize(e0[:,1]), normalize(e0[:,2]))

    # saved_positions=(true, true) because we discontinuously modify `u`.
    # This is independent of `saveat` and `save_everystep`
    cb = DiscreteCallback(scalingcondition, scale!; save_positions=(true, true))

    saved_values = SavedValues(Float64, ScaleRecord)
    scb = SavingCallback(savevalues, saved_values; save_everystep=false, saveat=zs, tdir=-1)

    p = WavefieldIntegrationParams(params.topheight, ea, frequency, bfield, species, params)

    # WARNING: Without `lazy=false` (necessary since we're using DiscreteCallback) don't
    # use continuous solution output! Also, we're only saving at zs.
    prob = ODEProblem{false}(dedz, e0, (first(zs), last(zs)), p)
    sol = solve(prob, solver; callback=CallbackSet(cb, scb),
                save_everystep=false, save_start=false, save_end=false,
                dt=dt, force_dtmin=force_dtmin, maxiters=maxiters,
                atol=tolerance, rtol=tolerance)

    if unscale
        e = unscalewavefields(saved_values)
    else
        records = saved_values.saveval
        e = [records[i].e for i in eachindex(records)]
    end

    return e
end

"""
    boundaryscalars(R, Rg, e1, e2, isotropic::Bool=false)
    boundaryscalars(R, Rg, e, isotropic::Bool=false)

Compute coefficients `(b1, b2)` required to sum the two wavefields vectors `e1` and `e2` or
both columns of `e` for the total wavefield at the ground as ``e = b1*e1 + b2*e2``.

The `b1` and `b2` that satisfy the waveguide boundary conditions are only valid for true
eigenangles of the waveguide.

!!! note

    This function assumes that the reflection coefficients `R` and `Rg` and the wavefield
    vectors `e1`, `e2` are at the ground.

    These boundary conditions
"""
function boundaryscalars(R, Rg, e1, e2, isotropic::Bool=false)
    # This is similar to excitation factor calculation, using the waveguide mode condition

    ex1, ey1, hx1, hy1 = e1[1], -e1[2], e1[3], e1[4]  # the ey component was originally -Ey
    ex2, ey2, hx2, hy2 = e2[1], -e2[2], e2[3], e2[4]

    # NOTE: Because `R`s and `e`s are at the ground, `hy0` and `ey0` are the same regardless
    # of `earthcurvature` being `true` or `false`.
    hy0 = 1  # == (1 + Rg[1,1])/(1 + Rg[1,1])
    ey0 = 1  # == (1 + Rg[2,2])/(1 + Rg[2,2])

    paral = 1 - R[1,1]*Rg[1,1]
    perp = 1 - R[2,2]*Rg[2,2]

    if isotropic
        if abs2(paral) < abs2(perp)
            b1 = hy0/hy1
            b2 = zero(hy1)
            # eyg = 0
        else
            b1 = zero(ey2)
            b2 = ey0/ey2
            # hyg = 0
        end
    else
        # polarization ratio Ey/Hy (often `f` or `fofr` in papers)
        if abs2(paral) < abs2(perp)
            EyHy = ((1 + Rg[2,2])*R[2,1]*Rg[1,1])/((1 + Rg[1,1])*perp)
        else
            EyHy = ((1 + Rg[2,2])*paral)/((1 + Rg[1,1])*R[1,2]*Rg[2,2])
        end
        a = (-ey1 + EyHy*hy1)/(ey2 - EyHy*hy2)
        hysum = hy1 + a*hy2
        b1 = hy0/hysum
        b2 = b1*a
    end

    return b1, b2
end

boundaryscalars(R, Rg, e, isotropic::Bool=false) =
    boundaryscalars(R, Rg, e[:,1], e[:,2], isotropic)

"""
    fieldstrengths!(EH, zs, me::ModeEquation; params=LMPParams())
    fieldstrengths!(EH, zs, ea::EigenAngle, frequency::Frequency, bfield::BField,
        species, ground::Ground; params=LMPParams())

Compute ``(Ex, Ey, Ez, Hx, Hy, Hz)ᵀ`` wavefields vectors as elements of `EH` by fullwave
integration at each height in `zs`.

The wavefields are scaled to satisfy the waveguide boundary conditions, which is only valid
at solutions of the mode equation.
"""
function fieldstrengths!(EH, zs, ea::EigenAngle, frequency::Frequency, bfield::BField,
    species, ground::Ground; params=LMPParams())

    @assert length(EH) == length(zs)
    zs[end] == 0 || @warn "Waveguide math assumes fields and reflection coefficients are
        calculated at the ground (`z = 0`)."

    e = integratewavefields(zs, ea, frequency, bfield, species; params=params)
    R = bookerreflection(ea, e[end])
    Rg = fresnelreflection(ea, ground, frequency)
    b1, b2 = boundaryscalars(R, Rg, e[end], isisotropic(bfield))

    S = ea.sinθ
    @inbounds for i in eachindex(EH)
        M = susceptibility(zs[i], frequency, bfield, species; params=params)

        # Scale to waveguide boundary conditions
        w = e[i][:,1]*b1 + e[i][:,2]*b2

        ex = w[1]
        ey = -w[2]
        ez = -(w[4]*S + M[3,1]*w[1] - M[3,2]*w[2])/(1 + M[3,3])
        hx = w[3]
        hy = w[4]
        hz = -w[2]*S

        EH[i] = SVector(ex, ey, ez, hx, hy, hz)
    end

    return nothing
end

function fieldstrengths!(EH, zs, me::ModeEquation; params=LMPParams())
    @unpack ea, frequency, waveguide = me
    @unpack bfield, species, ground = waveguide

    fieldstrengths!(EH, zs, ea, frequency, bfield, species, ground; params=params)
end

"""
    fieldstrengths(zs, me::ModeEquation; params=LMPParams())

Preallocate vector of wavefields `EH`, then call [`fieldstrengths!`](@ref) and return `EH`.

Each element of `EH` is an `SVector` of ``ex, ey, ez, hx, hy, hz``.
"""
function fieldstrengths(zs, me::ModeEquation; params=LMPParams())
    EH = Vector{SVector{6, ComplexF64}}(undef, length(zs))
    fieldstrengths!(EH, zs, me; params=params)
    return EH
end

"""
    calculate_wavefields!(wavefields, adjoint_wavefields, frequency, waveguide,
        adjoint_waveguide; params=LMPParams())

Compute fields of `wavefields` in-place scaled to satisfy the `waveguide` boundary
conditions.

This function implements the method of integrating wavefields suggested by
[Pitteway1965].

# References

[Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields, reflexion
    coefficients and polarizations for long radio waves in the lower ionosphere. I.,” Phil.
    Trans. R. Soc. Lond. A, vol. 257, no. 1079, pp. 219–241, Mar. 1965,
    doi: 10.1098/rsta.1965.0004.
"""
function calculate_wavefields!(wavefields, adjoint_wavefields, frequency, waveguide,
    adjoint_waveguide; params=LMPParams())

    @unpack bfield, species, ground = waveguide

    ground == adjoint_waveguide.ground ||
        throw(ArgumentError("waveguide and adjoint_waveguide should have same Ground"))
    species == adjoint_waveguide.species ||
        throw(ArgumentError("waveguide and adjoint_waveguide should have same Species"))

    adjoint_bfield = adjoint_waveguide.bfield

    zs = heights(wavefields)
    modes = eigenangles(wavefields)
    adjoint_modes = eigenangles(adjoint_wavefields)
    modes == adjoint_modes || @warn "Full mode conversion physics assumes adjoint
        wavefields are calculated for eigenangles of the original waveguide."

    for m in eachindex(modes)
        fieldstrengths!(view(wavefields,:,m), zs, modes[m], frequency, bfield, species,
                        ground; params=params)
        fieldstrengths!(view(adjoint_wavefields,:,m), zs, adjoint_modes[m], frequency,
                        adjoint_bfield, species, ground; params=params)
    end

    return nothing
end
