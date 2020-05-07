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

########
# Scenario

function scenario()
    bfield = BField(50e-6, deg2rad(68), deg2rad(111))
    tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
    ground = Ground(15, 0.001)

    electrons = Constituent(qₑ, mₑ,
                            z -> waitprofile(z, 75, 0.32),
                            electroncollisionfrequency)

    # electrons = Constituent(qₑ, mₑ,
    #         z -> waitprofile(z, 75, 0.32, H),
    #         z -> electroncollisionfrequency(z, H))

    # Resonant EigenAngle
    ea = EigenAngle(1.45964665843992 - 0.014974434753336im)

    ztop = LWMS.TOPHEIGHT
    zs = ztop:-100:zero(LWMS.TOPHEIGHT)

    return bfield, tx, ground, electrons, ea, zs
end


#==
modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)
tolerance = 1e-8

modes = LWMS.findmodes(origcoords, modeparams, tolerance)

ea = modes[argmax(real(modes))]  # largest real resonant mode
==#

########
# Does order of `q[1]`, `q[1]` matter in sharpboundaryreflection?
# Have to manually uncomment line in `_sharpboundaryreflection`

#==
tmp = LWMS.sharpboundaryreflection(ea, M)
tmpq = copy(LWMS.BOOKER_QUARTIC_ROOTS)

tmp2 = LWMS.sharpboundaryreflection(ea, M)
tmpq2 = copy(LWMS.BOOKER_QUARTIC_ROOTS)
==#

########
# Compare Booker quartic computed with M and T
# Runtime is dominated by `roots!`, but the version with `T` is slightly faster

bfield, tx, ground, electrons, ea, zs = scenario()

M = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
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

# Verify initialwavefields produces a valid solution
e = LWMS.initialwavefields(T)
q = LWMS.BOOKER_QUARTIC_ROOTS

for i = 1:2
    @test T*e[:,i] ≈ q[i]*e[:,i]
end

# Confirm `T` matrix vs dense array `T` both work
for i = 1:2
    @test T*e[:,i] ≈ Array(T)*e[:,i]
end


########
# Confirm reflection coefficient from wavefields at top height

bfield, tx, ground, electrons, ea, zs = scenario()

Mtop = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
Ttop = LWMS.tmatrix(ea, Mtop)
etop = LWMS.initialwavefields(Ttop)

wavefieldR = LWMS.vacuumreflectioncoeffs(ea, etop[:,1], etop[:,2])

# "modulus" (abs) of each component should be <=1
@test all(abs.(wavefieldR) .<= 1)
@test LWMS.sharpboundaryreflection(ea, Mtop) ≈ wavefieldR


########
# Integration with scaling

function integration_test()
    bfield, tx, ground, electrons, ea, zs = scenario()

    # Initial conditions
    Mtop = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
    Ttop = LWMS.tmatrix(ea, Mtop)
    e0 = LWMS.initialwavefields(Ttop)

    # Check normalized fields are still a valid solution
    e0n = hcat(normalize(e0[:,1]), normalize(e0[:,2]))

    q = LWMS.BOOKER_QUARTIC_ROOTS
    for i = 1:2
        @test Ttop*e0n[:,i] ≈ q[i]*e0n[:,i]
    end
    e0 = e0n

    cb = LWMS.DiscreteCallback(LWMS.scalingcondition, LWMS.scale!, save_positions=(true, true))
    saved_values = LWMS.SavedValues(eltype(zs), LWMS.ScaleRecord{eltype(zs), real(eltype(e0))})
    scb = LWMS.SavingCallback(LWMS.save_values, saved_values,
                              save_everystep=true, saveat=zs[2:end-1],
                              tdir=-1)

    p = LWMS.WavefieldIntegrationParams{eltype(e0)}(ea, tx.frequency, bfield, electrons)

    prob = ODEProblem{false}(LWMS.dedz, e0, (first(zs), last(zs)), p)
    sol = solve(prob, callback=CallbackSet(cb, scb),
                save_everystep=false, save_start=false, save_end=false,
                rtol=1e-8, atol=1e-8)


    return ea, saved_values
end

ea, saved_values = integration_test()
wavefieldRs = [LWMS.vacuumreflectioncoeffs(ea, s.e[:,1], s.e[:,2]) for s in saved_values.saveval]

# Compare to Budden integration of R
function dr_integration()
    bfield, tx, ground, electrons, ea, zs = scenario()

    modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)
    Mtop = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
    Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
    prob = ODEProblem{false}(LWMS.dRdz, Rtop, (first(zs), last(zs)), (ea, modeparams))
    sol = DifferentialEquations.solve(prob, Vern7(), abstol=1e-8, reltol=1e-8,
                                      saveat=saved_values.t, save_everystep=false)

    return sol
end

sol = dr_integration()

@test all(isapprox.(wavefieldRs, sol.u, atol=1e-7))


# Try again with unscaled fields
e = LWMS.unscalewavefields(saved_values)
wavefieldRs = [LWMS.vacuumreflectioncoeffs(ea, s[:,1], s[:,2]) for s in e]

@test all(isapprox.(wavefieldRs, sol.u, atol=1e-7))

function homogeneous_iono_test()
    # Check integration with homogeneous ionosphere
    # Integrate through homogeneous medium w/ sharp lower boundary and compare to
    # bookerquartic solution
    # See, e.g. Pitteway 1965 pg 234

    bfield, tx, ground, electrons, ea, zs = scenario()

    # LWMS.EARTHCURVATURE[] = false

    homogeneous = Constituent(qₑ, mₑ, z -> 1e8, z -> 5e6)

    homozs = 200e3:-500:50e3

    e, reczs = LWMS.integratewavefields(homozs, ea, tx.frequency, bfield, electrons)

    e1 = [s[:,1] for s in e]
    e2 = [s[:,2] for s in e]

    e1 = reshape(reinterpret(ComplexF64, e1), 4, :)
    e2 = reshape(reinterpret(ComplexF64, e2), 4, :)

    # plotly()
    # plot(real.(e1)', saved_values.t/1000)
    # plot(abs.(e2)', saved_values.t/1000)

    idx = findfirst(z->z==last(homozs), reczs)
    e1_idx = e1[:,idx]
    e2_idx = e2[:,idx]

    # need to scale to booker solution, which has second element == 1
    e1_idx *= 1/e1_idx[2]
    e2_idx *= 1/e2_idx[2]

    # Booker solution
    M = LWMS.susceptibility(last(homozs), tx.frequency, bfield, homogeneous)
    T = LWMS.tmatrix(ea, M)
    booker = LWMS.initialwavefields(T)

    q = LWMS.BOOKER_QUARTIC_ROOTS

    for i = 1:2
        @test T*booker[:,i] ≈ q[i]*booker[:,i]
    end

    # TODO: This test is supposed to work somehow...
    # plot(real(e2[1,:]), recz)

    @test_broken e1_idx ≈ booker[:,1]
    @test_broken e2_idx ≈ booker[:,2]

    # LWMS.EARTHCURVATURE[] = true
end

homogeneous_iono_test()


#==
Vacuum (ground) boundary condition
==#

e1, e2 = saved_values.saveval[end].e[:,1], saved_values.saveval[end].e[:,2]
R = vacuumreflectioncoeffs(ea, e1, e2)
Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)
b1, b2 = wavefieldboundary(R, Rg, e1, e2)


e = sol[:,end]
e1, e2 = e[:,1], e[:,2]

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


# `dt` argument is initial stepsize
# TODO: Use StepsizeLimiters in DifferentialEquations.jl to adaptively cap step
# size with λ in ionosphere / 20 or something like that




########

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
