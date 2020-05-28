using CSV
import DataFrames  # normalize conflicts with LinearAlgebra
using StaticArrays
using LinearAlgebra
using Interpolations

using DifferentialEquations
using Plots

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C


basepath = "/home/forrest/research/LAIR/ModeSolver/wavefields"

# Read in data from Piggott1965 plots
day_ne = CSV.read(joinpath(basepath,"piggot1965_day_ne.csv"),header=["ne","z"])
night_ne = CSV.read(joinpath(basepath,"piggot1965_night_ne.csv"),header=["ne","z"])
nu = CSV.read(joinpath(basepath,"piggot1965_nu.csv"),header=["nu","z"])


day_itp = interpolate((day_ne.z,), day_ne.ne, Gridded(Linear()))
day_expt = extrapolate(day_itp, Line())
dayne(z) = (tne = day_expt(z); tne > 0 ? tne : 0.001)


night_itp = interpolate((night_ne.z,), night_ne.ne, Gridded(Linear()))
night_expt = extrapolate(night_itp, Line())
nightne(z) = (tne = night_expt(z); tne > 0 ? tne : 0.001)


nu_itp = interpolate((nu.z,), nu.nu, Gridded(Linear()))
nu_expt = extrapolate(nu_itp, Line())
nuz(z) = (tne = nu_expt(z); tne > 0 ? tne : 0.001)

tzs = 0.0:100.0
plot(dayne.(tzs), tzs, xaxis=:log10)
plot(nightne.(tzs), tzs, xaxis=:log10)
plot(nuz.(tzs), tzs, xaxis=:log10)


########
# Pitteway 1965 fig 2 scenario

bfield = BField(50e-6, deg2rad(68), deg2rad(111))
tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(16e3), 100e3)
ground = Ground(15, 0.001)

electrons = Constituent(qₑ, mₑ, z->dayne(z/1000)*1e6, z->nuz(z/1000))

ea = EigenAngle(deg2rad(complex(40.0)))

ztop = 90e3
zs = ztop:-50:0.0

########
# Integrate wavefields

Mtop = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
Ttop = LWMS.tmatrix(ea, Mtop)
e0 = LWMS.initialwavefields(Ttop)

e0 = hcat(normalize(e0[:,1]), normalize(e0[:,2]))

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

e1 = [s.e[:,1] for s in saved_values.saveval]
e2 = [s.e[:,2] for s in saved_values.saveval]

e1 = reshape(reinterpret(ComplexF64, e1), 4, :)
e2 = reshape(reinterpret(ComplexF64, e2), 4, :)

escaled = LWMS.unscalewavefields(saved_values)

es1 = [t[:,1] for t in escaled]
es2 = [t[:,2] for t in escaled]

es1 = reshape(reinterpret(ComplexF64, es1), 4, :)
es2 = reshape(reinterpret(ComplexF64, es2), 4, :)

en, reczs = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, electrons)

en1 = [t[:,1] for t in en]
en2 = [t[:,2] for t in en]

en1 = reshape(reinterpret(ComplexF64, en1), 4, :)
en2 = reshape(reinterpret(ComplexF64, en2), 4, :)


df = DataFrames.DataFrame(z=saved_values.t/1000,
                          real_Ex1=real(es1[1,:]), imag_Ex1=imag(es1[1,:]), abs_Ex1=abs.(es1[1,:]),
                          real_Ey1=real(es1[2,:]), imag_Ey1=imag(es1[2,:]), abs_Ey1=abs.(es1[2,:]),
                          real_Ex2=real(es2[1,:]), imag_Ex2=imag(es2[1,:]), abs_Ex2=abs.(es2[1,:]),
                          real_Hx2=real(es2[3,:]), imag_Hx2=imag(es2[3,:]), abs_Hx2=abs.(es2[3,:]))

CSV.write(joinpath(basepath, "pitteway65_fig2.csv"), df)

plot(df.real_Ex1,df.z)
plot!(df.imag_Ex1,df.z)
plot!(df.abs_Ex1,df.z)
plot!(-df.abs_Ex1,df.z)
ylims!(50,80)

# Compare scaled and unscaled
df = DataFrames.DataFrame(z=saved_values.t/1000,
                          real_Ex1=real(es1[1,:]), real_Hx2=real(es2[3,:]),
                          real_Ex1u=real(e1[1,:]), real_Hx2u=real(e2[3,:]))

CSV.write(joinpath(basepath, "pitteway65_fig2_scaling.csv"), df)


# plot!(real(e1[1,:]),saved_values.t/1000)
# plot!(real(en1[1,:]),saved_values.t/1000)

#==
plotly()
# gr()

# e1
plot(abs.(e1)', saved_values.t/1000, color="red")
plot!(abs.(es1)', saved_values.t/1000, color="black")
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
==#


########
# Compare reflection coefficients

wavefieldRs = [LWMS.vacuumreflectioncoeffs(ea, s.e[:,1], s.e[:,2]) for s in saved_values.saveval]

modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)
Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
prob = ODEProblem{false}(LWMS.dRdz, Rtop, (first(zs), last(zs)), (ea, modeparams))
sol = DifferentialEquations.solve(prob, Vern7(), abstol=1e-8, reltol=1e-8,
                                  saveat=saved_values.t, save_everystep=false)


df = DataFrames.DataFrame(z=saved_values.t/1000,
               wf_R11=abs.(getindex.(wavefieldRs,1)),
               wf_R21=abs.(getindex.(wavefieldRs,2)),
               wf_R12=abs.(getindex.(wavefieldRs,3)),
               wf_R22=abs.(getindex.(wavefieldRs,4)),
               dr_R11=abs.(sol[1,:]),
               dr_R21=abs.(sol[2,:]),
               dr_R12=abs.(sol[3,:]),
               dr_R22=abs.(sol[4,:]))

CSV.write(joinpath(basepath, "pitteway65_fig2_refl.csv"), df)


########
# Pitteway 1965 fig 3 scenario

bfield = BField(50e-6, deg2rad(68), deg2rad(111))
tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(202e3), 100e3)
ground = Ground(15, 0.001)

# electrons = Constituent(qₑ, mₑ,
#                         z -> waitprofile(z, 80, 0.5),
#                         electroncollisionfrequency)

electrons = Constituent(qₑ, mₑ, z->nightne(z/1000)*1e6, z->nuz(z/1000))

ea = EigenAngle(deg2rad(complex(0.0)))

ztop = 110e3
zs = ztop:-10:0.0  # some very sharp transitions at this higher frequency

e, reczs = LWMS.integratewavefields(zs, ea, tx.frequency, bfield, electrons)

es1 = [t[:,1] for t in e]
es2 = [t[:,2] for t in e]

es1 = reshape(reinterpret(ComplexF64, es1), 4, :)
es2 = reshape(reinterpret(ComplexF64, es2), 4, :)

df = DataFrames.DataFrame(z=reczs/1000,
                          real_Ey1=real(es1[2,:]), imag_Ey1=imag(es1[2,:]), abs_Ey1=abs.(es1[2,:]),
                          real_Hx2=real(es2[3,:]), imag_Hx2=imag(es2[3,:]), abs_Hx2=abs.(es2[3,:]))

CSV.write(joinpath(basepath, "pitteway65_fig3.csv"), df)

plot(df.abs_Hx2,df.z)
plot!(-df.abs_Hx2,df.z)
plot!(df.real_Hx2,df.z)
plot!(df.imag_Hx2,df.z)
ylims!(74,102)


########
# Compare reflection coefficients

wavefieldRs = [LWMS.vacuumreflectioncoeffs(ea, t[:,1], t[:,2]) for t in e]

modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)
Mtop = LWMS.susceptibility(first(zs), tx.frequency, bfield, electrons)
Rtop = LWMS.sharpboundaryreflection(ea, Mtop)
prob = ODEProblem{false}(LWMS.dRdz, Rtop, (first(zs), last(zs)), (ea, modeparams))
sol = DifferentialEquations.solve(prob, Vern7(), abstol=1e-8, reltol=1e-8,
                                  saveat=reczs, save_everystep=false)


df = DataFrames.DataFrame(z=reczs/1000,
               wf_R11=abs.(getindex.(wavefieldRs,1)),
               wf_R21=abs.(getindex.(wavefieldRs,2)),
               wf_R12=abs.(getindex.(wavefieldRs,3)),
               wf_R22=abs.(getindex.(wavefieldRs,4)),
               dr_R11=abs.(sol[1,:]),
               dr_R21=abs.(sol[2,:]),
               dr_R12=abs.(sol[3,:]),
               dr_R22=abs.(sol[4,:]))

plotly()
plot(df.wf_R11,df.z)
plot!(df.dr_R11,df.z)

CSV.write(joinpath(basepath, "pitteway65_fig3_refl.csv"), df)
