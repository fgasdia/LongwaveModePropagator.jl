using Test
using LinearAlgebra
using StaticArrays
using ElasticArrays  # resizable multidimensional arrays

using OrdinaryDiffEq
using DifferentialEquations  # loading this to see what is chosen as default alg

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

tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
ground = Ground(15, 0.001)
bfield = BField(50e-6, deg2rad(50), deg2rad(300))
electrons = Constituent(qₑ, mₑ,
        z -> waitprofile(z, 75, 0.32, H),
        z -> electroncollisionfrequency(z, H))

# ea = modes[argmax(real(getfield.(modes, :θ)))]  # largest real resonant mode
ea = EigenAngle(complex(1.47152908, -0.0121998877))

# TODO: figure out how to confirm this ht is high enough
ztop = 120e3

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


########
#==
Where are the roots?

This shows that there are problems with magnetic dip angles of ±1°.

Also, the loop below is slow.
==#

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

########

#==
Integrate the wavefields
==#

function integratefields(e, p, z)
    M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)

    return LWMS.dedz(e, tx.frequency, T)
end


# First, regular integration
ea = EigenAngle(deg2rad(complex(40.0,0.0)))
bfield = BField(50e-6, deg2rad(68), deg2rad(111))
M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

# normalize e0
# normalizefields(e) = hcat(normalize(e[:,1]), normalize(e[:,2]))
# e0 = normalizefields(e0)

# `dt` argument is initial stepsize
prob = ODEProblem{false}(integratefields, e0, (ztop, 0.0))
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

#==
# Trying BigFloat, but this does no better
e0big = @SArray [parse(Complex{BigFloat}, string(e0[i,j])) for i in 1:4, j in 1:2]

prob = ODEProblem{false}(integratefields, e0big, BigFloat.((ztop, 0.0)))
solbig = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8)

@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'abs(e)'"
# @gp :- "set yrange [-10:10]"
for i = 1:4
    @gp :- real.(solbig[i,:]) solbig.t./1e3 "w l title '$(labels[i])'"
    @gp :- imag.(solbig[i,:]) solbig.t./1e3 "w l dt 2 title '$(labels[i])'"
end
==#

########
# Try rescaling
struct WavefieldScaling{T}
    # ortho::Array{T}
    e1norm::Vector{T}
    count::Vector{Int}
    # bnorm::Array{T}
end

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

"""
function scalewavefields(e1::AbstractVector, e2::AbstractVector, s::WavefieldScaling)
    # Orthogonalize vectors `e1` and `e2` (Gram-Schmidt process)
    # `dot` for complex vectors automatically conjugates first vector
    e1_dot_e1 = dot(e1, e1)  # == sum(abs2.(e1))
    a = dot(e1, e2)/e1_dot_e1
    e2 -= a*e1

    # Normalize `e1` and `e2`
    aterm = 1/sqrt(e1_dot_e1)
    e1 *= aterm
    e2 /= norm(e2)  # == dot(e2, e2)/sqrt(dot(e2, e2)) == normalize(e2)

    # p.ortho +=
    # s.anorm *= aterm
    # p.bnorm *=
    @inbounds s.count[1] += 1

    return e1, e2
end

"""
    myscale(e, s)

!!! note

    This function only applies scaling to the first 2 columns of `e`.
"""
function scalewavefields(e::AbstractArray, s::WavefieldScaling)
    e1, e2 = scalewavefields(e[:,1], e[:,2], s)

    # this works for a 4×2 `e` because `e[3:end]` will return a 4×0 array
    return hcat(e1, e2, e[:,3:end])
end

"""
    unscalewavefields(sol, p)
"""
function unscalewavefields(sol, s::WavefieldScaling)
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

function integratefieldswithscaling(e, p, z)
    # (u, p, t)

    z = p.z
    tx = p.tx
    bfield = p.bfield
    species = p.species
    s = p.s

    # Only want this at .t
    # Orthonomalize vectors
    e = scalewavefields(e, s)

    M = LWMS.susceptibility(z, tx.frequency, bfield, species)
    T = LWMS.tmatrix(ea, M)

    return LWMS.dedz(e, tx.frequency, T)
end


# Initial fields
ea = EigenAngle(complex(1.47152908, -0.0121998877))
M = LWMS.susceptibility(ztop, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

s = WavefieldScaling([],[0])
p = (z=ztop, tx=tx, bfield=bfield, species=electrons, s=s)

prob = ODEProblem(integratefieldswithscaling, e0, (ztop, 0.0), p)
sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8, dt=tx.frequency.λ/50)

zs = 0.0:1e3:ztop

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
