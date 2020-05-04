using Test
using LinearAlgebra
using StaticArrays
using ElasticArrays  # resizable multidimensional arrays
using Parameters

# using OrdinaryDiffEq
using DifferentialEquations  # loading this to see what is chosen as default alg

using Plots

using GRPF

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C

# To profile in Juno
# @profiler (for i = 1:1000; LWMS.fcn(); end)

#==
# Scenario
==#

bfield = BField(50e-6, deg2rad(68), deg2rad(111))
tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
ground = Ground(15, 0.001)

electrons = Constituent(qₑ, mₑ,
                        z -> waitprofile(z, 75, 0.32),
                        electroncollisionfrequency)

# electrons = Constituent(qₑ, mₑ,
#         z -> waitprofile(z, 75, 0.32, H),
#         z -> electroncollisionfrequency(z, H))

# ea = modes[argmax(real(getfield.(modes, :θ)))]  # largest real resonant mode
# ea = modes[argmax(real(modes))]
ea = EigenAngle(1.45964665843992 - 0.014974434753336im)

# TODO: figure out how to confirm this ht is high enough
ztop = 100e3

@with_kw struct WavefieldIntegrationParams{T1,T2,F,G}
    zprev::T1
    z::T1
    ortho_scalar::Complex{T2}
    e1_scalar::T2
    e2_scalar::T2
    frequency::Frequency
    bfield::BField
    species::Constituent{F,G}
end

#==
modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)
tolerance = 1e-8

modes = LWMS.findmodes(origcoords, modeparams, tolerance)
==#

########
# Does order of `q[1]`, `q[1]` matter in sharpboundaryreflection?
# Have to manually uncomment line in `_sharpboundaryreflection`

#==
ea = modes[argmax(real(getfield.(modes, :θ)))]
z = 120e3  # "TOP HT" TODO: figure out how to confirm this ht is high enough
M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

tmp = LWMS.sharpboundaryreflection(ea, M)
tmpq = copy(LWMS.BOOKER_QUARTIC_ROOTS)

tmp2 = LWMS.sharpboundaryreflection(ea, M)
tmpq2 = copy(LWMS.BOOKER_QUARTIC_ROOTS)
==#

########
# Compare Booker quartic computed with M and T
# NOTE: Runtime is dominated by `roots!` so currently no meaningful difference
# between the two

#==
M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)

LWMS.bookerquartic!(ea, M)
qM, BM = copy(LWMS.BOOKER_QUARTIC_ROOTS), copy(LWMS.BOOKER_QUARTIC_COEFFS)

LWMS.bookerquartic!(T)
qT, BT = copy(LWMS.BOOKER_QUARTIC_ROOTS), copy(LWMS.BOOKER_QUARTIC_COEFFS)

@test sort(qM, by=real) ≈ sort(qT, by=real)
@test sort(eigvals(Array(T)), by=real) ≈ sort(qT, by=real)  # eigvals is >20 times slower than bookerquartic

# Confirm Booker quartic is directly satisfied
for i in eachindex(qT)
    booker = qT[i]^4 + BT[4]*qT[i]^3 + BT[3]*qT[i]^2 + BT[2]*qT[i] + BT[1]
    @test isapprox(booker, 0, atol=1e-7)
end
==#

########
# We choose the 2 roots corresponding to upward travelling waves as being those that lie
# close to the positive real and negative imaginary axis (315° on the complex plane)

#==
e = LWMS.initialwavefields(T)
q = LWMS.BOOKER_QUARTIC_ROOTS

for i = 1:2
    @test T*e[:,i] ≈ q[i]*e[:,i]
end

for i = 1:2
    @test T*e[:,i] ≈ Array(T)*e[:,i]
end
==#

########
#==
Where are the roots?

This shows that there are problems with magnetic dip angles of ±1°.

Also, the loop below is slow.
==#

#==
# dips = -90:90
dips = vcat(-90:-2, 2:90) # Without ±1°
azims = 0:359

angles = [EigenAngle(complex(deg2rad(r),deg2rad(i))) for r = 30:89, i = -5:-1]
angles = vec(angles)

function calcroots(dips, azims, angles)
    M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
    qs = Array{MVector{4,ComplexF64}}(undef, length(angles), length(dips), length(azims))

    failindices = ElasticArray{Int}(undef, 3, 0)
    for az in eachindex(azims)
        for d in eachindex(dips)
            bfield = BField(50e-6, deg2rad(dips[d]), deg2rad(azims[az]))

            for a in eachindex(angles)
                LWMS.bookerquartic!(angles[a], M)
                q, B = copy(LWMS.BOOKER_QUARTIC_ROOTS), copy(LWMS.BOOKER_QUARTIC_COEFFS)

                for i in eachindex(q)
                    booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
                    if !isapprox(booker, 0, atol=1e-3)
                        append!(failindices, [a, d, az])
                    end
                end
                qs[a,d,az] = q
            end
        end
    end

    return qs
end

qs = calcroots(dips, azims, angles)

# unique(failindices[2,:])
# BUG: Worst performance is at dips[90] and dips[92] which corresponds to ±1°

qs1, qs2, qs3, qs4 = getindex.(qs,1), getindex.(qs,2), getindex.(qs,3), getindex.(qs,4)

# qr = hcat([vec(real(q)) for q in (qs1, qs2, qs3, qs4)]...)
# qi = hcat([vec(imag(q)) for q in (qs1, qs2, qs3, qs4)]...)

# We need to downsample because there's a large number of points for laptop
numels = reduce(*, size(qs))
randidxs = rand(1:numels, numels÷100)

#==
using Makie

scene = Scene();

maxx = 0
maxy = 0
for q in (qs1,)#, qs2, qs3, qs4)
    scatter!(scene, vec(real(q))[randidxs], vec(imag(q))[randidxs]);

    # Figure out axis lims
    qmaxx = max(abs(maximum(real(q))), abs(minimum(real(q))))
    qmaxy = max(abs(maximum(imag(q))), abs(minimum(imag(q))))

    global maxx
    global maxy
    qmaxx > maxx && (maxx = qmaxx)
    qmaxy > maxy && (maxy = qmaxy)
end
roundupto(n, x) = round(Int, (x+n/2)/n)*n
maxx = roundupto(5, maxx)
maxy = roundupto(5, maxy)

xlims!(scene, (-maxx, maxx));
ylims!(scene, (-maxy, maxy));

display(scene)
==#


# using GRUtils
# plot([vec(real(qs1))[randidxs] vec(real(qs2))[randidxs] vec(real(qs3))[randidxs] vec(real(qs4))[randidxs]],
#      [vec(imag(qs1))[randidxs] vec(imag(qs2))[randidxs] vec(imag(qs3))[randidxs] vec(imag(qs4))[randidxs]],
#      "x", markersize=0.5);
#
# xlabel("Re")
# ylabel("Im")
# savefig("roots.svg")

==#


########
# Integrate wavefields


zs = ztop:-100:0.0


#==
Calculate ionosphere reflection coefficient from `e`
==#

"""
    vacuumreflectioncoeffs(ea, e)

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
perpendicular to it, respectively [Budden1988] (pg 552).

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

See Budden 1988 ch 18 sec 7
"""
function vacuumreflectioncoeffs(ea::EigenAngle{T}, e1::AbstractArray{T2}, e2::AbstractArray{T2}) where {T,T2}
    C = ea.cosθ
    Cinv = 1/C

    # TODO: Special Sv matrix (also useful elsewhere?)
    Sv_inv = SMatrix{4,4,T,16}(Cinv, 0, -Cinv, 0,
                               0, -1, 0, -1,
                               0, -Cinv, 0, Cinv,
                               1, 0, 1, 0)

    f1 = Sv_inv*e1
    f2 = Sv_inv*e2

    out_type = promote_type(T, T2)
    U = SMatrix{2,2,out_type,4}(f1[1], f1[2], f2[1], f2[2])
    D = SMatrix{2,2,out_type,4}(f1[3], f1[4], f2[3], f2[4])

    return D/U
end

# Test R at "top"
Mtop = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
Ttop = LWMS.tmatrix(ea, Mtop)
etop = LWMS.initialwavefields(Ttop)
@test LWMS.sharpboundaryreflection(ea, Mtop) ≈ vacuumreflectioncoeffs(ea, etop[:,1], etop[:,2])

# BUG: broken wavefield calculation?
@test_broken all(abs.(vacuumR) .<= 1)  # Budden 88 says "modulus" (abs) of each component should be <=1

# Compare to Budden integration of R
modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)
Mtop = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
prob = ODEProblem{false}(LWMS.dRdz, Rtop, (ztop, 0.0), (ea, modeparams))
sol = DifferentialEquations.solve(prob, Vern7(), abstol=1e-8, reltol=1e-8)

# NOTE: this doesn't work with bigfloats, but it _does_ appear to work with scaling
@test_broken R ≈ sol[end]

#==
Vacuum (ground) boundary condition
==#

R = vacuumreflectioncoeffs(ea, e1, e2)
Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)

f = LWMS.modalequation(R, Rg)
@test isapprox(f, 0, atol=1e-3)  # is it worth iterating ea? dRdz is closer...

"""
    wavefieldboundary(R, Rg, e1, e2)

Calculate coefficients `b1`, `b2` required to sum `e1`, `e2` for total wavefield
as ``e = b1*e1 + b2*e2``.

This process is used by LWPC and is similar to the calculation of excitation
factor at the ground because it makes use of the mode condition.
"""
function wavefieldboundary(R, Rg, e1, e2)
    # This process is similar to excitation factor calculation, using the
    # waveguide mode condition

    Ex1, Ey1, Hx1, Hy1 = e1[1], -e1[2], e1[3], e1[4]
    Ex2, Ey2, Hx2, Hy2 = e2[1], -e2[2], e2[3], e2[4]

    # at the ground, modify for flat earth
    Hy0 = 1  # == (1 + Rg[1,1])/(1 + Rg[1,1])
    # ey0 = 1  # == (1 + Rg[2,2])/(1 + Rg[2,2])

    # polarization ratio Ey/Hy (often `f` or `fofr` in papers)
    abparal = 1 - Rg[1,1]*R[1,1]
    abperp = 1 - Rg[2,2]*R[2,2]
    if abs2(abparal) < abs2(abperp)
        pol = ((1 + Rg[2,2])*Rg[1,1]*R[2,1])/((1 + Rg[1,1])*(1 - Rg[2,2]*R[2,2]))
    else
        pol = ((1 + Rg[2,2])*abparal)/((1 + Rg[1,1])*Rg[2,2]*R[1,2])
    end

    a = (-Ey1 + pol*Hy1)/(Ey2 - pol*Hy2)

    Hysum = Hy1 + a*Hy2

    b1 = Hy0/Hysum
    b2 = b1*a

    return b1, b2
end

#==
Integration with bigfloats
==#

function integratefields(e, p, z)
    @unpack frequency, bfield, species = p

    M = LWMS.susceptibility(z, frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    return LWMS.dedz(e, frequency, T)
end


M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

# p = (;magnetoionic=(frequency=tx.frequency, bfield=bfield, species=electrons))
p = WavefieldIntegrationParams(0,0,complex(0),0,0,tx.frequency, bfield, electrons)
e0 = big.(e0)

# `dt` argument is initial stepsize
# TODO: Use StepsizeLimiters in DifferentialEquations.jl to adaptively cap step
# size with λ in ionosphere / 10 or something like that
# NOTE: The BigFloats make this quite slow...
prob = ODEProblem{false}(integratefields, e0, (ztop, 0.0), p)
sol = DifferentialEquations.solve(prob, abstol=1e-9, reltol=1e-9,
                                  dt=tx.frequency.λ/50, progress=true)

e = sol[:,end]
e1, e2 = e[:,1], e[:,2]

# need to make isbitstype for SArray
e1 = convert.(ComplexF64, e1)
e2 = convert.(ComplexF64, e2)

R = vacuumreflectioncoeffs(ea, e1, e2)
Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)
b1, b2 = wavefieldboundary(R, Rg, e1, e2)

totale = sol(zs)[:,1,:]*b1 + sol(zs)[:,2,:]*b2

plotly()
labels = ["Ex" "Ey" "Hx" "Hy"]
colors = palette(:tab10)[1:4]
p = plot(real.(totale'), zs/1000,
         label=labels, palette=palette(:tab10), linecolor=[1 2 3 4], linewidth=2)
plot!(p, imag.(totale'), zs/1000,
      ls=:dash, label=nothing, palette=palette(:tab10), linecolor=[1 2 3 4], linewidth=2)
xlabel!("e")
ylabel!("altitude (km)")


########
#==
Integration with scaling
==#

"""
    scalewavefields(e1, e2, s)

Orthonormalize vectors `e1` and `e2` and records scaling terms in `s`.

First applies Gram-Schmidt orthogonalization and then scales the vectors so they
each have length 1, i.e. `norm(e1) == norm(e2) == 1`. This is the technique
suggested by [^Pitteway1965] to counter numerical swamping during integration of
wavefields.

!!! note

    Elements of `s` are modified in place.


# References

[^Pitteway1965]: M. L. V. Pitteway, “The numerical calculation of wave-fields,
reflexion coefficients and polarizations for long radio waves in the lower
ionosphere. I.,” Phil. Trans. R. Soc. Lond. A, vol. 257, no. 1079,
pp. 219–241, Mar. 1965, doi: 10.1098/rsta.1965.0004.

[^Smith1974]: G. H. Smith and M. L. V. Pitteway, “Fortran Program for Obtaining
Wavefields of Penetrating, Nonpenetrating, and Whistler Modes of Radio Waves in
the Ionosphere,” in ELF-VLF Radio Wave Propagation, 1974, pp. 69–86.
"""
function scalewavefields(e1::AbstractVector, e2::AbstractVector)
    # Orthogonalize vectors `e1` and `e2` (Gram-Schmidt process)
    # `dot` for complex vectors automatically conjugates first vector
    e1_dot_e1 = real(dot(e1, e1))  # == sum(abs2.(e1)), NOTE: imag == 0
    a = dot(e1, e2)/e1_dot_e1  # purposefully unsigned
    e2 -= a*e1

    # Normalize `e1` and `e2`
    e1_scale_val = 1/sqrt(e1_dot_e1)
    e2_scale_val = 1/norm(e2)  # == 1/sqrt(dot(e2,e2))
    e1 *= e1_scale_val
    e2 *= e2_scale_val  # == normalize(e2)

    return e1, e2, a, e1_scale_val, e2_scale_val
end

"""
    scalewavefields(e, s)

!!! note

    This function only applies scaling to the first 2 columns of `e`.
"""
function scalewavefields(e::AbstractArray)
    e1, e2, a, e1_scale_val, e2_scale_val = scalewavefields(e[:,1], e[:,2])

    # this works for a 4×2 `e` because `e[3:end]` will return a 4×0 array
    e = hcat(e1, e2, e[:,3:end])

    return e, a, e1_scale_val, e2_scale_val
end

struct ScaleRecord{T1,T2}
    zprev::T1
    z::T1
    e::SMatrix{4,2,Complex{T2},8}
    ortho_scalar::Complex{T2}
    e1_scalar::T2
    e2_scalar::T2
end


# if any element of `e` has a real or imaginary component >= 1...
condition(e, z, integrator) = any(x -> (real(x) >= 1 || imag(x) >= 1), e)

function affect!(integrator)
    new_e, new_orthos, new_e1s, new_e2s = scalewavefields(integrator.u)

    #==
    NOTE: `integrator.t` is the "time" of the _proposed_ step. Therefore,
    integrator.t` might equal `0.0`, for example, before it's actually gotten
    to the bottom. Instead, `integrator.prevt` is the last `t` on the "left"
    side of the `integrator`, which covers the local interval [`tprev`, `t`].
    The "condition" is met at `integrator.t` and `integrator.t` is the time at
    which the affect occurs.
    However, it is not guaranteed that each tprev, t directly abuts the next tprev, t
    ==#

    println("affect! ($(integrator.tprev), $(integrator.t))")

    # Last set of scaling values
    @unpack frequency, bfield, species = integrator.p

    # NOTE: we must entirely reconstruct the entire NamedTuple from scratch
    integrator.p = WavefieldIntegrationParams(integrator.tprev, integrator.t,
                                              new_orthos*new_e1s/new_e2s,
                                              new_e1s, new_e2s,
                                              frequency, bfield, species)

    integrator.u = new_e

    # set_proposed_dt!(integrator, 0.2*get_proposed_dt(integrator))

    return nothing
end
# saved_positions=(true, true) because we discontinuously modify `u`. This is
# independent of saveat and save_everystep
cb = DiscreteCallback(condition, affect!, save_positions=(true, true))

saved_values = SavedValues(typeof(ztop), ScaleRecord{Float64, Float64})
save_values(u, t, integrator) = ScaleRecord(integrator.p.zprev,
                                            integrator.p.z,
                                            u,
                                            integrator.p.ortho_scalar,
                                            integrator.p.e1_scalar,
                                            integrator.p.e2_scalar)

# `save_everystep` because we need to make sure we save when affect! occurs
# `saveat=zs[2:end-1]` because otherwise we double save end points
scb = SavingCallback(save_values, saved_values,
                     save_everystep=true, saveat=zs[2:end-1],
                     tdir=-1)  # necessary because we're going down!

Mtop = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
Ttop = LWMS.tmatrix(ea, Mtop)
e0 = LWMS.initialwavefields(Ttop)
# TODO: Normalize e0?

# p = (z=ztop, ortho_scalar=zero(eltype(e0)), e1_scalar=1.0, e2_scalar=1.0,
     # magnetoionic=(frequency=tx.frequency, bfield=bfield, species=electrons))
p = WavefieldIntegrationParams(ztop, ztop, zero(eltype(e0)), 1.0, 1.0, tx.frequency, bfield, electrons)
prob = ODEProblem{false}(integratefields, e0, (ztop, zero(ztop)), p)
sol = solve(prob, callback=CallbackSet(cb, scb),
            save_everystep=false, save_start=false, save_end=false)

# plot(real(sol[:,1,:])', sol.t/1000)

e1, e2 = saved_values.saveval[end].e[:,1], saved_values.saveval[end].e[:,2]
R = vacuumreflectioncoeffs(ea, e1, e2)
Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)
b1, b2 = wavefieldboundary(R, Rg, e1, e2)

#==
# Compare to Budden integration of R
modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)
Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
prob = ODEProblem{false}(LWMS.dRdz, Rtop, (ztop, 0.0), (ea, modeparams))
sol = DifferentialEquations.solve(prob, Vern7(), abstol=1e-8, reltol=1e-8)

@test R ≈ sol[end]
==#

"""
    unscalewavefields(sol, p)

I think it works like this:

- Scaling is applied to each level from the top down
- The bottom level does not get unscaled. Instead, we will be referencing
the above levels to this level.
- The level above the bottom level needs to be additionally scaled by the amount
that was applied to get to the bottom level.
- The next level up (2 above the bottom level) needs to be scaled by the amount
applied to the next level and then the bottom level, i.e. we need to keep track
of a cumulative correction on the way back up.

We do not need to keep track of cumulative values on the way down.
We only need to update the cumulative correction on the way up when there is a
new correction.
"""
function unscalewavefields!(e::AbstractVector, saved_values::SavedValues)

    z = saved_values.t
    records = saved_values.saveval

    osum = zero(records[1].ortho_scalar)
    prod_e1 = one(records[1].e1_scalar)
    prod_e2 = one(records[1].e2_scalar)

    # ref_zprev = first(records).zprev
    # ref_z = first(records).z  # purposefully first(records) != last(records)
    ref_zprev = last(records).zprev
    ref_z = last(records).z

    # Initialize ref_z
    prod_e1 *= last(records).e1_scalar
    prod_e2 *= last(records).e2_scalar

    # Unscaling we go from the bottom up
    for i in reverse(eachindex(e))
        # Unpack variables
        record_zprev = records[i].zprev
        record_z = records[i].z
        scaled_e = records[i].e
        ortho_scalar = records[i].ortho_scalar
        e1_scalar = records[i].e1_scalar
        e2_scalar = records[i].e2_scalar

        # Only update ref_z when there is both:
        # 1) we have reached the height where a new ref_z should go into effect
        # 2) there is a new ref_z
        if z[i] > record_z && record_z > ref_z
            ref_z = record_z

            prod_e1 *= e1_scalar
            prod_e2 *= e2_scalar
        end

        if i == lastindex(e)  # == [end]
            # Bottom doesn't require correction
            e[i] = scaled_e
        else
            # e2 = (scaled_e[i][:,2] - osum*scaled_e[i][:,1])*prod_e2
            # e2 = (scaled_e[:,2] - ortho_cumulative_scalar[i]*scaled_e[:,1])*e2_cumulative_scalar[i]
            e2 = scaled_e[:,2]*prod_e2
            e1 = scaled_e[:,1]*prod_e1
            e[i] = hcat(e1,e2)
        end

        # # if zs[i] >= ref_zprev
        # if record_z != ref_z
        #     # osum *= e1_cumulative_scalar[i]/e2_cumulative_scalar[i]
        #     # osum += ortho_cumulative_scalar[i]
        #     prod_e1 *= e1_scalar
        #     prod_e2 *= e2_scalar
        #
        #     ref_zprev = record_zprev
        #     ref_z = record_z
        # end
    end

    return nothing
end

function unscalewavefields(saved_values::SavedValues)
    # Array of SArray{Tuple{4,2}, Complex{Float64}}
    e = Vector{typeof(saved_values.saveval[1].e)}(undef, length(saved_values.saveval))

    unscalewavefields!(e, saved_values)

    return e
end

e = unscalewavefields(saved_values)

e1 = [s.e[:,1] for s in saved_values.saveval]
e2 = [s.e[:,2] for s in saved_values.saveval]

e1 = reshape(reinterpret(ComplexF64, e1), 4, :)
e2 = reshape(reinterpret(ComplexF64, e2), 4, :)

ne1 = [t[:,1] for t in e]
ne2 = [t[:,2] for t in e]

ne1 = reshape(reinterpret(ComplexF64, ne1), 4, :)
ne2 = reshape(reinterpret(ComplexF64, ne2), 4, :)

affect_zs = [s.z for s in saved_values.saveval]

plotly()
# gr()
plot(abs.(e1)', saved_values.t/1000, color="red")
plot!(abs.(ne1)', saved_values.t/1000, color="black")
scatter!(zeros(length(sol.t)),sol.t/1000, markersize=6, markercolor=nothing, markerstrokecolor="blue")
scatter!(zeros(length(affect_zs)), affect_zs/1000,
        markershape=:rect, markersize=3, markercolor=nothing, markerstrokecolor="black")
vline!([1])




function fieldstrengths()

end


########

#==
Integrate the wavefields - TEST
==#

function integratefields(e, p, z)
    frequency, bfield, species = p

    M = LWMS.susceptibility(z, frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    return LWMS.dedz(e, frequency, T)
end


# First, regular integration
ea = EigenAngle(deg2rad(complex(40.0,0.0)))
bfield = BField(50e-6, deg2rad(68), deg2rad(111))
M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)
q = copy(LWMS.BOOKER_QUARTIC_ROOTS)  # generated in `initialwavefields`


# normalize e0
# This seemed like a good idea, but the wavefields still blow up
#==
normalizefields(e) = hcat(normalize(e[:,1]), normalize(e[:,2]))
e0 = normalizefields(e0)

for i = 1:2
    @test T*e0[:,i] ≈ q[i]*e0[:,i]
end
==#


p = (frequency=tx.frequency, bfield=bfield, species=electrons)

# `dt` argument is initial stepsize
# TODO: Use StepsizeLimiters in DifferentialEquations.jl to adaptively cap step
# size with λ in ionosphere / 20 or something like that
prob = ODEProblem{false}(integratefields, e0, (ztop, 0.0), p)
sol = solve(prob, abstol=1e-8, reltol=1e-8, dt=tx.frequency.λ/50)

# `solve` chose composite Vern9 and Rodas5
# sol.alg

zs = 0.0:100:120e3

using Gnuplot
@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'e'"
# @gp :- "set yrange [-10:10]"

labels = Dict(1=>"Ex", 2=>"Ey", 3=>"Hx", 4=>"Hy")
for i = 1:4
    @gp :- real.(sol(zs)[i,:]) zs./1e3 "w l title '$(labels[i])'"
    @gp :- imag.(sol(zs)[i,:]) zs./1e3 "w l dt 2 title '$(labels[i])'"
end

# Clearly, they blew up

#
# Trying BigFloat, but this does no better
e0big = @SArray [parse(Complex{BigFloat}, string(e0[i,j])) for i in 1:4, j in 1:2]

prob = ODEProblem{false}(integratefields, e0big, BigFloat.((ztop, 0.0)), p)
solbig = solve(prob, abstol=1e-18, reltol=1e-18)

@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'e₁'"
# @gp :- "set yrange [-10:10]"
for i = 1:4
    @gp :- real.(solbig[i,:]) solbig.t./1e3 "w l lc '$(colors[i])' title '$(labels[i])'"
    @gp :- imag.(solbig[i,:]) solbig.t./1e3 "w l dt 2 lc '$(colors[i])' title '$(labels[i])'"
end


@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'e₂'"
# @gp :- "set yrange [-10:10]"
for i = 5:8
    @gp :- real.(solbig[i,:]) solbig.t./1e3 "w l lc '$(colors[i-4])' title '$(labels[i-4])'"
    @gp :- imag.(solbig[i,:]) solbig.t./1e3 "w l dt 2 lc '$(colors[i-4])' title '$(labels[i-4])'"
end
==#

########
# Try rescaling

mutable struct WavefieldMatrix{T,T2} <: DEDataMatrix{T}
    # Extra fields are carried through solver, and can be accessed as e.g. `u.count`

    # TEMP XXX This should be an SMatrix
    # x::SMatrix{4,2,T,8}  # must be called `x` (I think...)
    x::Matrix{T}
    e1norm::Vector{T2}
    count::Int
end


"""
    unscalewavefields(sol, p)
"""
function unscalewavefields(sol)
    # Back substitution of normalizing values applied during integration

    # Array of SArray{Tuple{4,2}, Complex{Float64}}
    e = Array{typeof(first(sol))}(undef, length(sol))

    osum = 0
    proda = 1
    prodb = 1

    for i in eachindex(sol)
        # e[:,2,i] = prodb*(sol[:,2,1] - osum*sol[:,1,1])
        e[:,1,i] = proda*sol[:,1,1]

        # osum = osum*anorm/bnorm + ortho
        proda *= p.anorm
        # prodb *= bnorm
    end

    # For performance, should be something like this
    #==
    function normSoA(x)
      out = 0.0
      @simd for i in 1:length(x.real)
        @inbounds out += sqrt(x.real[i]^2 + x.imag[i]^2)
      end
      out
    end
    ==#

end

function scalingcondition(e, z, integrator)
    # When `scalingcondition()` == 0, trigger `scalewavefields!()`

    norm(e[:,1]) > 1000 || abs2(dot(e[:,1], e[:,2])) > 1 ? true : false
end



function integratefields!(de, e, p, z)
    tx, bfield, species = p

    M = LWMS.susceptibility(z, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    de .= LWMS.dedz(e, tx.frequency, T)
end


# Initial fields
ea = EigenAngle(complex(1.47152908, -0.0121998877))
M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

p = (tx=tx, bfield=bfield, species=electrons)

# TEMP
# u0 = WavefieldMatrix(e0, Vector{ComplexF64}(), 0)
u0 = WavefieldMatrix(Array(e0), Vector{ComplexF64}(), 0)

e0big = [parse(Complex{BigFloat}, string(e0[i,j])) for i in 1:4, j in 1:2]
u0big = WavefieldMatrix(e0big, Vector{Complex{BigFloat}}(), 0)

# FunctionCallingCallback doesn't only evaluate at sol.t... it evaluates every fcn call
# tdir = -1 b/c t[1] > t[end]
# fcc = FunctionCallingCallback((e, z, integrator)->scalewavefields(e, s), func_everystep=true, tdir=-1)


function affect!(integrator)
    # Need to loop through `full_cache` to ensure that all internal caches are
    # also updated.
    for c in full_cache(integrator)
        x, e1norm = scalewavefields(c.x)
        c.x .= x
        c.e1norm = e1norm
        c.count += 1
    end

    # TODO: Also look at updating dt when e is rescaled

    return nothing
end

cb = DiscreteCallback(scalingcondition, affect!, save_positions=(true,true))

prob = ODEProblem(integratefields!, e0big, (ztop, 0.0), p)
solbig = solve(prob, callback=cb, abstol=1e-18, reltol=1e-18, dt=tx.frequency.λ/50)


prob = ODEProblem(integratefields!, u0, (ztop, 0.0), p)
sol = solve(prob, callback=cb, abstol=1e-8, reltol=1e-8, dt=tx.frequency.λ/50)

# TODO: Use a callback on independence of e1 and e2 and then save at the points
# where callback is applied


zs = 0.0:250:ztop

using Gnuplot
@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'e'"
# @gp :- "set yrange [-10:10]"
labels = Dict(1=>"Ex", 2=>"Ey", 3=>"Hx", 4=>"Hy")
for i = 1:4
    @gp :- real.(sol[i,:]) sol.t./1e3 "w l title '$(labels[i])'"
    @gp :- imag.(sol[i,:]) sol.t./1e3 "w l dt 2 title '$(labels[i])'"
end

########
# Nagano et al 1975 method

# BUG?: Nagano et al 1975 specifies that abs(q[1]) > abs(q[2]), BUT my sorting
# method doesn't agree with this.

"""
z -> ztop:0.0

Technically, we don't need to compute the `e`s as we go. If we save the `K`
matrices, then
``e₁(zⱼ) = ∏ₛ₌ⱼ₊₁ᵐ⁻¹ (Kₛ) e₁(zₘ₋₁)``
``e₂₀(zⱼ) = aⱼ e₁(zⱼ) + ∏ₛ₌ⱼ₊₁ᵐ⁻¹ (Kₛ) e₂(zₘ₋₁)``

This might be useful if we only need values at specific heights.
"""
function nagano_integration(z, ea, frequency, bfield, species)
    if z[end] > z[1]
        @warn "nagano_integration should proceed downwards. z[end] $(z[end]) > z[1] $(z[1])"
    end

    e1 = Vector{SVector{4,ComplexF64}}(undef, length(z))
    e2 = Vector{SVector{4,ComplexF64}}(undef, length(z))
    avec = Vector{ComplexF64}(undef, length(z))
    Bvec = Vector{SMatrix{4,4,ComplexF64,16}}(undef, length(z))

    # `...cscalar` means _cumulative_ scalar values
    e1cscalar = Vector{ComplexF64}(undef, length(z))
    e2cscalar = Vector{ComplexF64}(undef, length(z))
    orthocscalar = Vector{ComplexF64}(undef, length(z))
    # TODO: check these types

    # Initialize cumulative vectors
    e1cscalar[1] = 1
    e2cscalar[1] = 1
    orthocscalar[1] = 0

    # cumulative_e1_scalars = one(ComplexF64)
    # cumulative_e2_scalar = one(ComplexF64)
    # cumulative_ortho_scalar = zero(ComplexF64)

    k = frequency.k

    # Temporary mutable matrix
    B = MMatrix{4,4,ComplexF64,16}(undef)
    B[2,:] = 1  # TODO: Special BMatrix type

    K = MMatrix{4,4,ComplexF64,16}(undef)

    for j in eachindex(z)

        if j > firstindex(z)
            tmp_e1 = K*e1[j-1]  # this is `K(j-1)` b/c it's from last loop
            tmp_e2 = K*e2[j-1]

            # Scale wavefields
            e1s, e2s, a, e1_scalar, e2_scalar = scalewavefields(tmp_e1, tmp_e2)

            e1[j] = e1s  # scaled `e1`
            e2[j] = e2s  # scaled `e2`
            avec[j] = a  # unsigned `a`
            # e1scale[j] = e1_scalar

            # First update scaled `a`, then update scale values
            orthocscalar[j] = orthocscalar[j-1] + a*(e1cscalar[j-1]/e2cscalar[j-1])
            e1cscalar[j] = e1cscalar[j-1]*e1_scalar
            e2cscalar[j] = e2cscalar[j-1]*e2_scalar
        end

        # Calculate `Kⱼ` for use in next step, remember e(zⱼ₋₁) = Kⱼ e(zⱼ)
        M = LWMS.susceptibility(z[j], frequency, bfield, species)
        T = LWMS.tmatrix(ea, M)
        LWMS.bookerquartic!(T)
        LWMS.sortquarticroots!(LWMS.BOOKER_QUARTIC_ROOTS)
        q = LWMS.BOOKER_QUARTIC_ROOTS

        for i = 1:4
            # TODO: Is this equivalent to LWMS.initialwavefields?
            # Apparently not? A little suspicious...
            # BUG: Try with my solution? (initialwavefields)
            α = (T[4,2]*T[3,2] - (T[3,2] - q[i]^2)*(T[4,4] - q[i]))/
                (T[3,1]*(T[4,4] - q[i]) - T[4,1]*T[3,4])
            β = (T[1,2] + α*(T[1,1] - q[i]))/T[1,4]

            B[1,i] = α
            # B[2,i] = 1
            B[3,i] = q[i]
            B[4,i] = β
        end

        if j == firstindex(z)
            # Initialize wavefields (at top height "m")
            e1[j] = SVector(B[:,1])  # == B*[1; 0; 0; 0]
            e1cscalar[j] = 1  # no scaling applied
            e2[j] = SVector(B[:,2])  # == B*[0; 1; 0; 0]
        end

        if j == lastindex(z)
            nothing
        else
            zstep = z[j+1] - z[j]  # next lower level - current level
            ξ = q*zstep
            Δ = SDiagonal(exp.(-im*k*ξ))
            K .= B*Δ/B
        end

        # Store values
        Bvec[j] = SMatrix(B)
    end

    #==
    # Reconstruct the wavefields from the bottom up
    # `aprod` is prod from s=1 (bottom height) to height j-1
    e1prod = one(eltype(e1cscalar))
    e2prod = one(eltype(e2cscalar))
    orthosum = zero(eltype(orthocscalar))
    for j in reverse(eachindex(z))
        if j == lastindex(z)  # == first(reverse(eachindex(z)))
            # don't unscale at bottom height
            continue
        end
        # e1prod *= avec[j+1]  # at previous (lower) height
        e1prod *= e1cscalar[j]
        e2prod *= e2cscalar[j]
        orthosum = orthosum*e1cscalar[j]/e2cscalar[j] + orthocscalar[j]
        e2[j] = (e2[j] - orthosum*e1[j])*e2prod
        e1[j] *= e1prod
        # e2[j] += aprod*e1[end]  # `e1` at bottom height
        # e1[j] *= e1scale[j]  # unnormalize `e1`, which was scaled after `e2`
    end
    ==#

    return e1, e2
end

# @test_skip vacuum wavefields  # (Nagano BC)


# BUG: Handle when Ne = 0
# electrons = Constituent(qₑ, mₑ,
#         z -> waitprofile(z, 75, 0.32, H),
#         z -> electroncollisionfrequency(z, H))
electrons = Constituent(qₑ, mₑ,
        z -> waitprofile(z, 75, 0.32),
        electroncollisionfrequency)

ea = EigenAngle(deg2rad(complex(40.0,0.0)))
bfield = BField(50e-6, deg2rad(68), deg2rad(111))
tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
ground = Ground(15, 0.001)

#==
electrons = Constituent(qₑ, mₑ,
        z -> waitprofile(z, 75, 0.32),
        electroncollisionfrequency)

ea = EigenAngle(deg2rad(complex(60.0,0.0)))
bfield = BField(50e-6, deg2rad(59), deg2rad(90))
tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(10e3), 100e3)
==#

zs = ztop:-100:0.0
e1, e2 = nagano_integration(zs, ea, tx.frequency, bfield, electrons)

e = e1 + e2


using Gnuplot
@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'e₁'"
# @gp :- "set yrange [-10:10]"
labels = Dict(1=>"Ex", 2=>"Ey", 3=>"Hx", 4=>"Hy")
colors = Dict(1=>"red", 2=>"blue", 3=>"green", 4=>"black")
for i = 1:4
    @gp :- real.(getindex.(e1,i)) zs/1e3 "w l lc '$(colors[i])' title '$(labels[i])'"
    @gp :- imag.(getindex.(e1,i)) zs/1e3 "w l dt 2 lc '$(colors[i])' title '$(labels[i])'"
end

@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'e₂'"
# @gp :- "set yrange [-10:10]"
labels = Dict(1=>"Ex", 2=>"Ey", 3=>"Hx", 4=>"Hy")
colors = Dict(1=>"red", 2=>"blue", 3=>"green", 4=>"black")
for i = 1:4
    @gp :- real.(getindex.(e2,i)) zs/1e3 "w l lc '$(colors[i])' title '$(labels[i])'"
    @gp :- imag.(getindex.(e2,i)) zs/1e3 "w l dt 2 lc '$(colors[i])' title '$(labels[i])'"
end


@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'e'"
# @gp :- "set yrange [-10:10]"
labels = Dict(1=>"Ex", 2=>"Ey", 3=>"Hx", 4=>"Hy")
for i = 1:4
    @gp :- real.(getindex.(e,i)) zs/1e3 "w l title '$(labels[i])'"
    @gp :- imag.(getindex.(e,i)) zs/1e3 "w l dt 2 title '$(labels[i])'"
end



########
# Confirm Nagano BookerQuartic solution is valid

M = LWMS.susceptibility(72e3, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
LWMS.bookerquartic!(T)
LWMS.sortquarticroots!(LWMS.BOOKER_QUARTIC_ROOTS)
q = LWMS.BOOKER_QUARTIC_ROOTS

B = MMatrix{4,4,ComplexF64,16}(undef)
B[2,:] = 1  # TODO: Special BMatrix type
for i = 1:4
    α = (T[4,2]*T[3,2] - (T[3,2] - q[i]^2)*(T[4,4] - q[i]))/
        (T[3,1]*(T[4,4] - q[i]) - T[4,1]*T[3,4])
    β = (T[1,2] + α*(T[1,1] - q[i]))/T[1,4]

    B[1,i] = α
    # B[2,i] = 1
    B[3,i] = q[i]
    B[4,i] = β
end


for i = 1:4
    @test_broken T*B[:,i] ≈ q[i]*B[:,i]
end

o = LWMS.initialwavefields(T)


########
#==
TESTS

and LWPC comparisons
==#

function homogeneousiono()
    # Accuracy check for integration of wavefields. See Pitteway1965 pg 234
    # Integrate through homogeneous medium with sharp lower boundary and compare
    # to bookerquartic solution
end
@test_skip homogeneousiono() ≈ bookerquartic()




function lwpcreflectioncoeffs(ea::EigenAngle, e1, e2)
    # From wf_r_mtrx.for

    C = ea.cosθ

    g12 = e1[1]*e2[2] - e2[1]*e1[2]
    g13 = e1[1]*e2[3] - e2[1]*e1[3]
    g14 = e1[1]*e2[4] - e2[1]*e1[4]
    g23 = e1[2]*e2[3] - e2[2]*e1[3]
    g24 = e1[2]*e2[4] - e2[2]*e1[4]
    g34 = e1[3]*e2[4] - e2[3]*e1[4]

    den = -g13 + C*(g34 - g12 + C*g24)

    d11 = g13 + C*(g34 + g12 + C*g24)
    d22 = g13 + C*(-g34 - g12 + C*g24)
    d12 = 2*C*g14
    d21 = 2*C*g23

    return SMatrix{2,2,eltype(den),4}(d11/den, d21/den, d12/den, d22/den)
end

vacuumR = vacuumreflectioncoeffs(ea, e1[end], e2[end])
@test vacuumR ≈ lwpcreflectioncoeffs(ea, e1[end], e2[end])



function lwpcscale(p1, p2)
    e1, e2 = MVector(p1), MVector(p2)

    # aterm → dot(e1, e1)
    aterm = 0
    for i = 1:4
        aterm += abs2(e1[i])
    end

    # term → dot(e1, e2)
    term = 0
    for i = 1:4
        term += conj(e1[i])*e2[i]
    end

    # term → dot(e1, e2)/dot(e1, e1)
    term /= aterm

    # e2 → e2 - dot(e1, e2)/dot(e1, e1)
    for i = 1:4
        e2[i] -= term*e1[i]
    end

    # bterm → dot(e2, e2)
    bterm = 0
    for i = 1:4
        bterm += abs2(e2[i])
    end

    # Normalize both vectors
    aterm = 1/sqrt(aterm)
    bterm = 1/sqrt(bterm)
    for i = 1:4
        e1[i] *= aterm
        e2[i] *= bterm
    end

    return SVector(e1), SVector(e2)
end

########
