using Test
using LinearAlgebra
using StaticArrays

using OrdinaryDiffEq

using GRPF

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C

# To profile in Juno
# @profiler (for i = 1:1000; LWMS.fcn(); end)

# Scenario

tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
ground = Ground(15, 0.001)
bfield = BField(50e-6, deg2rad(50), deg2rad(300))
electrons = Constituent(qₑ, mₑ,
        z -> waitprofile(z, 75, 0.32, H),
        z -> electroncollisionfrequency(z, H))

modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)
tolerance = 1e-6

modes = LWMS.findmodes(origcoords, modeparams, tolerance)

# Does order of `q[1]`, `q[1]` matter in sharpboundaryreflection?
# Have to manually uncomment line in `_sharpboundaryreflection`
ea = modes[argmax(real(getfield.(modes, :θ)))]
z = 120e3  # "TOP HT" TODO: figure out how to confirm this ht is high enough
M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

tmp = LWMS.sharpboundaryreflection(ea, M)
tmpq = copy(LWMS.BOOKER_QUARTIC_ROOTS)

tmp2 = LWMS.sharpboundaryreflection(ea, M)
tmpq2 = copy(LWMS.BOOKER_QUARTIC_ROOTS)

########
# Compare Booker quartic computed with M and T
# NOTE: Runtime is dominated by `roots!` so currently no meaningful difference
# between the two

# ea = modes[argmax(real(getfield.(modes, :θ)))]
ea = EigenAngle(complex(1.47152908, -0.0121998877))
z = 120e3  # "TOP HT" TODO: figure out how to confirm this ht is high enough
M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
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
==#

z = 120e3  # "TOP HT"
# dips = -90:90
dips = vcat(-90:-2, 2:90) # Without ±1°
azims = 0:359

angles = [EigenAngle(complex(deg2rad(r),deg2rad(i))) for r = 30:89, i = -5:-1]
angles = vec(angles)

qs = Array{MVector{4,ComplexF64}}(undef, length(angles), length(dips), length(azims))

using ElasticArrays
failindices = ElasticArray{Int}(undef, 3, 0)

for az in eachindex(azims)
    for d in eachindex(dips)
        bfield = BField(50e-6, deg2rad(dips[d]), deg2rad(azims[az]))

        for a in eachindex(angles)
            M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)

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
# unique(failindices[2,:])
# BUG: Worst performance is at dips[90] and dips[92] which corresponds to ±1°

using Gnuplot

qs1, qs2, qs3, qs4 = getindex.(qs,1), getindex.(qs,2), getindex.(qs,3), getindex.(qs,4)

numels = reduce(*, size(qs))
randidxs = rand(1:numels, numels÷50)

@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
@gp :- "unset key"
# @gp :- "set yrange [-10:10]"
for q in (qs1, qs2, qs3, qs4)
    @gp :- vec(real(q))[randidxs] vec(imag(q))[randidxs] "w p pt 7 ps 0.5"
end
@gp

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

# First, regular integration

function integratefields(e, p, z)
    M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)

    return LWMS.dedz(e, tx.frequency, T)
end

# Initial fields
ea = EigenAngle(complex(1.47152908, -0.0121998877))
z = 120e3  # "TOP HT" TODO: figure out how to confirm this ht is high enough
M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

prob = ODEProblem{false}(integratefields, e0, (120e3, 0.0))
sol = solve(prob, BS5(), abstol=1e-6, reltol=1e-6)

zs = 0.0:1e3:120e3

using Gnuplot
@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'abs(e)'"
# @gp :- "set yrange [-10:10]"
for i = 1:4
    @gp :- abs.(sol[i,:]) sol.t./1e3 "w l title '$i'"
end

# Clearly, they blew up

# Use big floats
ea = EigenAngle(complex(1.47152908, -0.0121998877))
z = 120e3  # "TOP HT" TODO: figure out how to confirm this ht is high enough
M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

e0big = @SArray [parse(Complex{BigFloat}, string(e0[i,j])) for i in 1:4, j in 1:2]

prob = ODEProblem{false}(integratefields, e0big, BigFloat.((120e3, 0.0)))
sol2 = solve(prob, BS5(), abstol=1e-6, reltol=1e-6)

@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'abs(e)'"
# @gp :- "set yrange [-10:10]"
for i = 1:4
    @gp :- abs.(sol2[i,:]) sol2.t./1e3 "w l title '$i'"
end

# This does not fix the problem?

# Try rescaling
function integratefieldswithscaling(e, p, z)
    # if something
        # `dot` for complex vectors automatically conjugates first vector
        a = -dot(e[:,1], e[:,2])/dot(e[:,1], e[:,1])
    # end
    e = hcat(e[:,1], e[:,2]+a*e[:,1])

    M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
    T = LWMS.tmatrix(ea, M)

    return LWMS.dedz(e, tx.frequency, T)
end


# Initial fields
ea = EigenAngle(complex(1.47152908, -0.0121998877))
z = 120e3  # "TOP HT" TODO: figure out how to confirm this ht is high enough
M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
e0 = LWMS.initialwavefields(T)

prob = ODEProblem(integratefieldswithscaling, e0, (120e3, 0.0))
sol = solve(prob, BS5(), abstol=1e-6, reltol=1e-6)

zs = 0.0:1e3:120e3

using Gnuplot
@gp "set auto fix"
@gp :- "set grid"
@gp :- "set offsets graph .05, graph .05, graph .05, graph .05"
# @gp :- "unset key"
@gp :- "set ylabel 'height (km)'" "set xlabel 'abs(e)'"
# @gp :- "set yrange [-10:10]"
for i = 1:4
    @gp :- abs.(sol[i,:]) sol.t./1e3 "w l title '$i'"
end
