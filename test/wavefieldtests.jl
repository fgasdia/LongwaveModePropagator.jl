using Test
using LinearAlgebra
using StaticArrays
using ElasticArrays  # resizable multidimensional arrays
using Parameters

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

ztop = 100e3


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
# Runtime is dominated by `roots!`, but the version with `T` is slightly faster

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


########
# We choose the 2 roots corresponding to upward travelling waves as being those that lie
# close to the positive real and negative imaginary axis (315° on the complex plane)

e = LWMS.initialwavefields(T)
q = LWMS.BOOKER_QUARTIC_ROOTS

for i = 1:2
    @test T*e[:,i] ≈ q[i]*e[:,i]
end

for i = 1:2
    # Confirm `T` matrix vs dense array `T` both work
    @test T*e[:,i] ≈ Array(T)*e[:,i]
end


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

M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

p = LWMS.WavefieldIntegrationParams(0, complex(0), 0, 0, ea, tx.frequency, bfield, electrons)

# `dt` argument is initial stepsize
# TODO: Use StepsizeLimiters in DifferentialEquations.jl to adaptively cap step
# size with λ in ionosphere / 10 or something like that
prob = ODEProblem{false}(integratefields, e0, (ztop, 0.0), p)
sol = DifferentialEquations.solve(prob, abstol=1e-8, reltol=1e-8,
                                  dt=tx.frequency.λ/50, progress=true)

e = sol[:,end]
e1, e2 = e[:,1], e[:,2]


########
#==
Integration with scaling
==#


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

# e1
plot(abs.(e1)', saved_values.t/1000, color="red")
plot!(abs.(ne1)', saved_values.t/1000, color="black")
scatter!(zeros(length(sol.t)),sol.t/1000, markersize=6, markercolor=nothing, markerstrokecolor="blue")
scatter!(zeros(length(affect_zs)), affect_zs/1000,
        markershape=:rect, markersize=3, markercolor=nothing, markerstrokecolor="black")
vline!([1])

# e2
plot(abs.(e2)', saved_values.t/1000, color="red")
plot!(abs.(ne2)', saved_values.t/1000, color="black")
scatter!(zeros(length(sol.t)),sol.t/1000, markersize=6, markercolor=nothing, markerstrokecolor="blue")
scatter!(zeros(length(affect_zs)), affect_zs/1000,
        markershape=:rect, markersize=3, markercolor=nothing, markerstrokecolor="black")
vline!([1])




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

R = vacuumreflectioncoeffs(ea, e1, e2)
Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)
b1, b2 = wavefieldboundary(R, Rg, e1, e2)



function fieldstrengths()

end


########

#==
Integrate the wavefields - TEST
==#


# `dt` argument is initial stepsize
# TODO: Use StepsizeLimiters in DifferentialEquations.jl to adaptively cap step
# size with λ in ionosphere / 20 or something like that
prob = ODEProblem{false}(integratefields, e0, (ztop, 0.0), p)
sol = solve(prob, abstol=1e-8, reltol=1e-8, dt=tx.frequency.λ/50)


########

# @test_skip vacuum wavefields  # (Nagano BC)


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
