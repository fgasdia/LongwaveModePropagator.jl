using LinearAlgebra
using StaticArrays
using DataFrames
using Gadfly
using Parameters
using Printf

# include("..\\src\\LongwaveModeSolver.jl")
using LongwaveModeSolver
const LWMS = LongwaveModeSolver

finitediff(v, h) = diff(v)./h
symmetricdiff(v, h) = (v[3:end]-v[1:end-2])./(2*h)
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

#==
Bookerquartic and initial R
==#
ω = 2π*24e3
z₀ = 50e3
z = 75e3
spec = LWMS.Constituent(-LWMS.fundamentalcharge, LWMS.electronmass,
    h -> LWMS.waitprofile(h, 75, 0.3), LWMS.electroncollisionfrequency)
bfield = LWMS.BField(0.5e-4, -0.2, -0.3, -0.9)
M = LWMS.susceptibility(ω, z₀, z, spec, bfield)

h = 0.5*π/180
eas = [LWMS.EigenAngle(th) for th in complex.(range(0.0, π/2; step=h), 10π/180)]
qBs = [LWMS.bookerquartic(ea, M) for ea in eas]
qs = getindex.(qBs, 1)
Bs = getindex.(qBs, 2)
sort!.(qs, by=LWMS.anglefrom315)

Ss = getfield.(eas, :sinθ)
Cs = getfield.(eas, :cosθ)
C²s = getfield.(eas, :cos²θ)

θs = getfield.(eas, :θ)

xticks = deg2rad.(range(0.0, 100.0; step=15))

# B
function dB(S, C, C², B1, M)
    dS = C
    dC² = -2*S*C
    dB3 = dS*(M[1,3] + M[3,1])
    # dB2 = dC²*(-M[3,3] - 1) - dC²*(M[1,1] + 1)
    dB2 = -dC²*(2 + M[1,1] + M[3,3])
    # dB1 = dC²*(-M[1,3] - M[3,1])*S + dS*(M[1,2]*M[2,3] + M[2,1]*M[3,2] -
        # (M[1,3] + M[3,1])*(M[1,2] + C²))
    dB1 = dS/S*B1 - S*dC²*(M[1,3] + M[3,1])
    # dB0 = dC²*(-M[1,2]*M[2,1] - M[1,3]*M[3,1] + (M[1,1] + 1)*(M[2,2] + C²) +
    #     (M[1,1] + 1) * (M[3,3] + C²))
    dB0 = dC²*(2*C²*(1 + M[1,1]) + M[3,3] + M[2,2] + M[1,1]*(M[3,3] + M[2,2]) -
        M[1,3]*M[3,1] - M[1,2]*M[2,1])

    return (dB0, dB1, dB2, dB3)
end

dBs = [dB(Ss[i], Cs[i], C²s[i], Bs[i][2], M) for i in 1:length(Ss)]

realBs = [real.(b) for b in Bs]
realdBs = [real.(db) for db in dBs]
imagBs = [imag.(b) for b in Bs]
imagdBs = [imag.(db) for db in dBs]

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:5
    tmpdf[!,Symbol("B$(i-1)")] = getfield.(Tuple.(realBs), i)
end
realdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=5), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:5
    tmpdf[!,Symbol("B$(i-1)")] = getfield.(Tuple.(imagBs), i)
end
imagdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=5), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dB$(i-1)")] = getfield.(realdBs, i)
    # fdrealB = finitediff(getfield.(Tuple.(realBs), i), deg2rad(1.0))
    # tmpdf[!,Symbol("finite_dB$(i-1)")] = [fdrealB; fdrealB[end]]
end
tmpdf[!,Symbol("dB4")] = fill(0.0, length(θs))
realddf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=5), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dB$(i-1)")] = getfield.(imagdBs, i)
    # fdimagB = finitediff(getfield.(Tuple.(imagBs), i), deg2rad(1.0))
    # tmpdf[!,Symbol("finite_dB$(i-1)")] = [fdimagB; fdimagB[end]]
end
tmpdf[!,Symbol("dB4")] = fill(0.0, length(θs))
imagddf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=5), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])


pr = plot(realdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realddf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagddf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\Bs.svg", 12inch, 7inch)

#########
# q
dq(B, dB, q) = -(((dB[4]*q + dB[3])*q + dB[2])*q + dB[1]) / (((4*B[5]*q + 3*B[4])*q + 2*B[3])*q + B[2])

dq_1 = dq.(Bs, dBs, getfield.(Tuple.(qs), 1))
dq_2 = dq.(Bs, dBs, getfield.(Tuple.(qs), 2))

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("q1")] = real.(getfield.(Tuple.(qs), 1))
tmpdf[!,Symbol("q2")] = real.(getfield.(Tuple.(qs), 2))

realqdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:q1, :q2])[!,1],
    val=stack(tmpdf[!,2:end], [:q1, :q2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("q1")] = imag.(getfield.(Tuple.(qs), 1))
tmpdf[!,Symbol("q2")] = imag.(getfield.(Tuple.(qs), 2))

imagqdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:q1, :q2])[!,1],
    val=stack(tmpdf[!,2:end], [:q1, :q2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dq1")] = real.(dq_1)
tmpdf[!,Symbol("dq2")] = real.(dq_2)
fdq_1 = finitediff(real.(getfield.(Tuple.(qs), 1)), h)
fdq_2 = finitediff(real.(getfield.(Tuple.(qs), 2)), h)
tmpdf[!,Symbol("finite_dq_1")] = [fdq_1; fdq_1[end]]
tmpdf[!,Symbol("finite_dq_2")] = [fdq_2; fdq_2[end]]

realdqdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dq1")] = imag.(dq_1)
tmpdf[!,Symbol("dq2")] = imag.(dq_2)
fdq_1 = finitediff(imag.(getfield.(Tuple.(qs), 1)), h)
fdq_2 = finitediff(imag.(getfield.(Tuple.(qs), 2)), h)
tmpdf[!,Symbol("finite_dq_1")] = [fdq_1; fdq_1[end]]
tmpdf[!,Symbol("finite_dq_2")] = [fdq_2; fdq_2[end]]

imagdqdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realqdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagqdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdqdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdqdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\qs.svg", 12inch, 7inch)

########
# Auxiliaries (Δ, T, P)

function dPTdelta(ea, M)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ

    @unpack q, B, D12, D32, D33, D11_1, D13_1, D31_1, Δ_1, invΔ_1, P_1, T_1,
        D11_2, D13_2, D31_2, Δ_2, invΔ_2, P_2, T_2, Δ, invΔ = LWMS._common_sharplyboundedR(ea, M)

    # Additional calculations required for dR/dθ
    dS = C
    dC = -S
    dC² = -2*S*C
    dB3 = dS*(M[1,3] + M[3,1])
    dB2 = -dC²*(2 + M[1,1] + M[3,3])
    dB1 = dS/S*B[2] - S*dC²*(M[1,3] + M[3,1])
    dB0 = dC²*(2*C²*(1 + M[1,1]) + M[3,3] + M[2,2] + M[1,1]*(M[3,3] + M[2,2]) -
        M[1,3]*M[3,1] - M[1,2]*M[2,1])

    dq_1 = -(((dB3*q[1] + dB2)*q[1] + dB1)*q[1] + dB0) /
        (((4*B[5]*q[1] + 3*B[4])*q[1] + 2*B[3])*q[1] + B[2])
    dq_2 = -(((dB3*q[2] + dB2)*q[2] + dB1)*q[2] + dB0) /
        (((4*B[5]*q[2] + 3*B[4])*q[2] + 2*B[3])*q[2] + B[2])

    dD33 = dC²

    dD11_1 = -2*q[1]*dq_1
    dD13_1 = dq_1*S + q[1]*dS
    dD31_1 = dD13_1  # dq_1*S + q[1]*dS

    dΔ_1 = dD11_1*D33 + D11_1*dD33 - dD13_1*D31_1 - D13_1*dD31_1
    dinvΔ_1 = -dΔ_1/Δ_1^2

    dP_1 = (-D12*dD33 + dD13_1*D32)*invΔ_1 + (D13_1*D32 - D12*D33)*dinvΔ_1
    dT_1 = dq_1*P_1 + q[1]*dP_1 -
        dS*(-D11_1*D32 + D12*D31_1)*invΔ_1 -
        S*(-dD11_1*D32 + D12*dD31_1)*invΔ_1 -
        S*(-D11_1*D32 + D12*D31_1)*dinvΔ_1

    dD11_2 = -2*q[2]*dq_2
    dD13_2 = dq_2*S + q[2]*dS
    dD31_2 = dD13_2  # dq_2*S + q[2]*dS

    dΔ_2 = dD11_2*D33 + D11_2*dD33 - dD13_2*D31_2 - D13_2*dD31_2
    dinvΔ_2 = -dΔ_2/Δ_2^2

    dP_2 = (-D12*dD33 + dD13_2*D32)*invΔ_2 + (D13_2*D32 - D12*D33)*dinvΔ_2
    dT_2 = dq_2*P_2 + q[2]*dP_2 -
        dS*(-D11_2*D32 + D12*D31_2)*invΔ_2 -
        S*(-dD11_2*D32 + D12*dD31_2)*invΔ_2 -
        S*(-D11_2*D32 + D12*D31_2)*dinvΔ_2

    dΔ = dT_1*C² + T_1*dC² + dT_1*C*q[2] + T_1*dC*q[2] + T_1*C*dq_2 + dP_1*C + P_1*dC +
        dP_1*q[2] + P_1*dq_2 - (dT_2*C² + T_2*dC²) -
        (dT_2*C*q[1] + T_2*dC*q[1] + T_2*C*dq_1) -
        (dP_2*C + P_2*dC) - (dP_2*q[1] + P_2*dq_1)
    dinvΔ = -dΔ/Δ^2

    return (Δ_1, dΔ_1, Δ_2, dΔ_2, P_1, dP_1, P_2, dP_2, invΔ_1, dinvΔ_1, invΔ_2, dinvΔ_2,
        D11_1, dD11_1, D11_2, dD11_2, Δ, dΔ, invΔ, dinvΔ, T_1, dT_1, T_2, dT_2)
end

auxs = [dPTdelta(ea, M) for ea in eas]

realΔ_1, imagΔ_1 = unzip(reim.(getfield.(auxs, 1)))
realdΔ_1, imagdΔ_1 = unzip(reim.(getfield.(auxs, 2)))
realΔ_2, imagΔ_2 = unzip(reim.(getfield.(auxs, 3)))
realdΔ_2, imagdΔ_2 = unzip(reim.(getfield.(auxs, 4)))
realP_1, imagP_1 = unzip(reim.(getfield.(auxs, 5)))
realdP_1, imagdP_1 = unzip(reim.(getfield.(auxs, 6)))
realP_2, imagP_2 = unzip(reim.(getfield.(auxs, 7)))
realdP_2, imagdP_2 = unzip(reim.(getfield.(auxs, 8)))
realinvΔ_1, imaginvΔ_1 = unzip(reim.(getfield.(auxs, 9)))
realdinvΔ_1, imagdinvΔ_1 = unzip(reim.(getfield.(auxs, 10)))
realinvΔ_2, imaginvΔ_2 = unzip(reim.(getfield.(auxs, 11)))
realdinvΔ_2, imagdinvΔ_2 = unzip(reim.(getfield.(auxs, 12)))
realD11_1, imagD11_1 = unzip(reim.(getfield.(auxs, 13)))
realdD11_1, imagdD11_1 = unzip(reim.(getfield.(auxs, 14)))
realD11_2, imagD11_2 = unzip(reim.(getfield.(auxs, 15)))
realdD11_2, imagdD11_2 = unzip(reim.(getfield.(auxs, 16)))
realΔ, imagΔ = unzip(reim.(getfield.(auxs, 17)))
realdΔ, imagdΔ = unzip(reim.(getfield.(auxs, 18)))
realinvΔ, imaginvΔ = unzip(reim.(getfield.(auxs, 19)))
realdinvΔ, imagdinvΔ = unzip(reim.(getfield.(auxs, 20)))
realT_1, imagT_1 = unzip(reim.(getfield.(auxs, 21)))
realdT_1, imagdT_1 = unzip(reim.(getfield.(auxs, 22)))
realT_2, imagT_2 = unzip(reim.(getfield.(auxs, 23)))
realdT_2, imagdT_2 = unzip(reim.(getfield.(auxs, 24)))

# D11_1 and D11_2
tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("D11_1")] = realD11_1
tmpdf[!,Symbol("D11_2")] = realD11_2

realD11df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:D11_1, :D11_2])[!,1],
    val=stack(tmpdf[!,2:end], [:D11_1, :D11_2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("D11_1")] = imagD11_1
tmpdf[!,Symbol("D11_2")] = imagD11_2

imagD11df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:D11_1, :D11_2])[!,1],
    val=stack(tmpdf[!,2:end], [:D11_1, :D11_2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dD11_1")] = realdD11_1
tmpdf[!,Symbol("dD11_2")] = realdD11_2
fdD11_1 = finitediff(realD11_1, h)
fdD11_2 = finitediff(realD11_2, h)
tmpdf[!,Symbol("finite_dD11_1")] = [fdD11_1; fdD11_1[end]]
tmpdf[!,Symbol("finite_dD11_2")] = [fdD11_2; fdD11_2[end]]

realdD11df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dD11_1")] = imagdD11_1
tmpdf[!,Symbol("dD11_2")] = imagdD11_2
fdD11_1 = finitediff(imagD11_1, h)
fdD11_2 = finitediff(imagD11_2, h)
tmpdf[!,Symbol("finite_dD11_1")] = [fdD11_1; fdD11_1[end]]
tmpdf[!,Symbol("finite_dD11_2")] = [fdD11_2; fdD11_2[end]]

imagdD11df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realD11df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagD11df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdD11df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdD11df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\D11s.svg", 12inch, 7inch)

# Δ_1 and Δ_2
tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("Δ_1")] = realΔ_1
tmpdf[!,Symbol("Δ_2")] = realΔ_2
tmpdf[!,Symbol("Δ")] = realΔ

realΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=3), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("Δ_1")] = imagΔ_1
tmpdf[!,Symbol("Δ_2")] = imagΔ_2
tmpdf[!,Symbol("Δ")] = imagΔ

imagΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=3), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dΔ_1")] = realdΔ_1
tmpdf[!,Symbol("dΔ_2")] = realdΔ_2
tmpdf[!,Symbol("dΔ")] = realdΔ

realdΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=3), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dΔ_1")] = imagdΔ_1
tmpdf[!,Symbol("dΔ_2")] = imagdΔ_2
tmpdf[!,Symbol("dΔ")] = imagdΔ

imagdΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=3), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\Δs.svg", 12inch, 7inch)

# invΔ_1 and invΔ_2
tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("invΔ_1")] = realinvΔ_1
tmpdf[!,Symbol("invΔ_2")] = realinvΔ_2
tmpdf[!,Symbol("invΔ")] = realinvΔ

realinvΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=3), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("invΔ_1")] = imaginvΔ_1
tmpdf[!,Symbol("invΔ_2")] = imaginvΔ_2
tmpdf[!,Symbol("invΔ")] = imaginvΔ

imaginvΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=3), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dinvΔ_1")] = realdinvΔ_1
tmpdf[!,Symbol("dinvΔ_2")] = realdinvΔ_2
tmpdf[!,Symbol("dinvΔ")] = realdinvΔ
fdinvΔ_1 = finitediff(realinvΔ_1, h)
fdinvΔ_2 = finitediff(realinvΔ_2, h)
fdinvΔ = finitediff(realinvΔ, h)
tmpdf[!,Symbol("finite_dinvΔ_1")] = [fdinvΔ_1; fdinvΔ_1[end]]
tmpdf[!,Symbol("finite_dinvΔ_2")] = [fdinvΔ_2; fdinvΔ_2[end]]
tmpdf[!,Symbol("finite_dinvΔ")] = [fdinvΔ; fdinvΔ[end]]

realdinvΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=6), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dinvΔ_1")] = imagdinvΔ_1
tmpdf[!,Symbol("dinvΔ_2")] = imagdinvΔ_2
tmpdf[!,Symbol("dinvΔ")] = imagdinvΔ
fdinvΔ_1 = finitediff(imaginvΔ_1, h)
fdinvΔ_2 = finitediff(imaginvΔ_2, h)
fdinvΔ = finitediff(imaginvΔ, h)
tmpdf[!,Symbol("finite_dinvΔ_1")] = [fdinvΔ_1; fdinvΔ_1[end]]
tmpdf[!,Symbol("finite_dinvΔ_2")] = [fdinvΔ_2; fdinvΔ_2[end]]
tmpdf[!,Symbol("finite_dinvΔ")] = [fdinvΔ; fdinvΔ[end]]

imagdinvΔdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=6), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realinvΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imaginvΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdinvΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdinvΔdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\invΔs.svg", 12inch, 7inch)

# P_1 and P_2
tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("P_1")] = realP_1
tmpdf[!,Symbol("P_2")] = realP_2

realPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:P_1, :P_2])[!,1],
    val=stack(tmpdf[!,2:end], [:P_1, :P_2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("P_1")] = imagP_1
tmpdf[!,Symbol("P_2")] = imagP_2

imagPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:P_1, :P_2])[!,1],
    val=stack(tmpdf[!,2:end], [:P_1, :P_2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dP_1")] = realdP_1
tmpdf[!,Symbol("dP_2")] = realdP_2
fdP_1 = finitediff(realP_1, h)
fdP_2 = finitediff(realP_2, h)
tmpdf[!,Symbol("finite_dP_1")] = [fdP_1; fdP_1[end]]
tmpdf[!,Symbol("finite_dP_2")] = [fdP_2; fdP_2[end]]

realdPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dP_1")] = imagdP_1
tmpdf[!,Symbol("dP_2")] = imagdP_2
fdP_1 = finitediff(imagP_1, h)
fdP_2 = finitediff(imagP_2, h)
tmpdf[!,Symbol("finite_dP_1")] = [fdP_1; fdP_1[end]]
tmpdf[!,Symbol("finite_dP_2")] = [fdP_2; fdP_2[end]]

imagdPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\Ps.svg", 12inch, 7inch)

# T_1 and T_2
tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("T_1")] = realT_1
tmpdf[!,Symbol("T_2")] = realT_2

realPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:T_1, :T_2])[!,1],
    val=stack(tmpdf[!,2:end], [:T_1, :T_2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("T_1")] = imagT_1
tmpdf[!,Symbol("T_2")] = imagT_2

imagPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=2), var=stack(tmpdf[!,2:end], [:T_1, :T_2])[!,1],
    val=stack(tmpdf[!,2:end], [:T_1, :T_2])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dT_1")] = realdT_1
tmpdf[!,Symbol("dT_2")] = realdT_2
fdT_1 = finitediff(realT_1, h)
fdT_2 = finitediff(realT_2, h)
tmpdf[!,Symbol("finite_dT_1")] = [fdT_1; fdT_1[end]]
tmpdf[!,Symbol("finite_dT_2")] = [fdT_2; fdT_2[end]]

realdPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
tmpdf[!,Symbol("dT_1")] = imagdT_1
tmpdf[!,Symbol("dT_2")] = imagdT_2
fdT_1 = finitediff(imagT_1, h)
fdT_2 = finitediff(imagT_2, h)
tmpdf[!,Symbol("finite_dT_1")] = [fdT_1; fdT_1[end]]
tmpdf[!,Symbol("finite_dT_2")] = [fdT_2; fdT_2[end]]

imagdPdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdPdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (deg)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\Ts.svg", 12inch, 7inch)

########
# Rs
RdRs = [LWMS.sharplybounded_R_dRdθ(ea, M) for ea in eas]
Rs = getfield.(RdRs, 1)
dRs = getfield.(RdRs, 2)

realRs = Tuple.(real.(Rs))
realdRs = Tuple.(real.(dRs))

imagRs = Tuple.(imag.(Rs))
imagdRs = Tuple.(imag.(dRs))

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("R$i")] = getfield.(realRs, i)
end
realRdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("R$i")] = getfield.(imagRs, i)
end
imagRdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dR$i")] = getfield.(realdRs, i)
end
for i = 1:4
    realfdR = finitediff(getfield.(realRs, i), h)
    tmpdf[!,Symbol("finite_dR$i")] = [realfdR; realfdR[end]]
    # realfdR = symmetricdiff(getfield.(realRs, i), h)
    # tmpdf[!,Symbol("finite_dR$i")] = [realfdR[1]; realfdR; realfdR[end]]
end

realdRdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dR$i")] = getfield.(imagdRs, i)
end
for i = 1:4
    imagfdR = finitediff(getfield.(imagRs, i), h)
    tmpdf[!,Symbol("finite_dR$i")] = [imagfdR; imagfdR[end]]
    # imagfdR = symmetricdiff(getfield.(imagRs, i), h)
    # tmpdf[!,Symbol("finite_dR$i")] = [imagfdR[1]; imagfdR; imagfdR[end]]
end

imagdRdf = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realRdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagRdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdRdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdRdf, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\Rs.svg", 12inch, 7inch)


#===
S Matrix
===#

Ts = [LWMS.tmatrix(ea, M) for ea in eas]
Ss = [LWMS.smatrix(ea, T) for (ea, T) in zip(eas, Ts)]
dSs = [LWMS.dsmatrixdθ(ea, M, T) for (ea, T) in zip(eas, Ts)]

S1s = getfield.(Ss, 1)
S2s = getfield.(Ss, 2)
S3s = getfield.(Ss, 3)
S4s = getfield.(Ss, 4)

dS1s = getfield.(dSs, 1)
dS2s = getfield.(dSs, 2)
dS3s = getfield.(dSs, 3)
dS4s = getfield.(dSs, 4)

realS1s, imagS1s = unzip(reim.(S1s))
realS2s, imagS2s = unzip(reim.(S2s))
realS3s, imagS3s = unzip(reim.(S3s))
realS4s, imagS4s = unzip(reim.(S4s))

realdS1s, imagdS1s = unzip(reim.(dS1s))
realdS2s, imagdS2s = unzip(reim.(dS2s))
realdS3s, imagdS3s = unzip(reim.(dS3s))
realdS4s, imagdS4s = unzip(reim.(dS4s))

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S11_$i")] = getfield.(Tuple.(realS1s), i)
end
realS1df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S11_$i")] = getfield.(Tuple.(imagS1s), i)
end
imagS1df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS11_$i")] = getfield.(Tuple.(realdS1s), i)
end
for i = 1:4
    realfdR = finitediff(getfield.(Tuple.(realS1s), i), h)
    tmpdf[!,Symbol("finite_dS11_$i")] = [realfdR; realfdR[end]]
end
realdS1df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS11_$i")] = getfield.(Tuple.(imagdS1s), i)
end
for i = 1:4
    imagfdR = finitediff(getfield.(Tuple.(imagS1s), i), h)
    tmpdf[!,Symbol("finite_dS11_$i")] = [imagfdR; imagfdR[end]]
end
imagdS1df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realS1df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagS1df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdS1df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdS1df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\S1s.svg", 12inch, 7inch)


# S12
tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S12_$i")] = getfield.(Tuple.(realS2s), i)
end
realS2df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S12_$i")] = getfield.(Tuple.(imagS2s), i)
end
imagS2df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS12_$i")] = getfield.(Tuple.(realdS2s), i)
end
for i = 1:4
    realfdR = finitediff(getfield.(Tuple.(realS2s), i), h)
    tmpdf[!,Symbol("finite_dS12_$i")] = [realfdR; realfdR[end]]
end
realdS2df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS12_$i")] = getfield.(Tuple.(imagdS2s), i)
end
for i = 1:4
    imagfdR = finitediff(getfield.(Tuple.(imagS2s), i), h)
    tmpdf[!,Symbol("finite_dS12_$i")] = [imagfdR; imagfdR[end]]
end
imagdS2df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realS2df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagS2df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdS2df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdS2df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\S2s.svg", 12inch, 7inch)


# S21
tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S21_$i")] = getfield.(Tuple.(realS3s), i)
end
realS3df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S21_$i")] = getfield.(Tuple.(imagS3s), i)
end
imagS3df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS21_$i")] = getfield.(Tuple.(realdS3s), i)
end
for i = 1:4
    realfdR = finitediff(getfield.(Tuple.(realS3s), i), h)
    tmpdf[!,Symbol("finite_dS21_$i")] = [realfdR; realfdR[end]]
end
realdS3df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS21_$i")] = getfield.(Tuple.(imagdS3s), i)
end
for i = 1:4
    imagfdR = finitediff(getfield.(Tuple.(imagS3s), i), h)
    tmpdf[!,Symbol("finite_dS21_$i")] = [imagfdR; imagfdR[end]]
end
imagdS3df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realS3df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagS3df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdS3df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdS3df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\S3s.svg", 12inch, 7inch)

# S22
tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S22_$i")] = getfield.(Tuple.(realS4s), i)
end
realS4df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("S22_$i")] = getfield.(Tuple.(imagS4s), i)
end
imagS4df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=4), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS22_$i")] = getfield.(Tuple.(realdS4s), i)
end
for i = 1:4
    realfdR = finitediff(getfield.(Tuple.(realS4s), i), h)
    tmpdf[!,Symbol("finite_dS22_$i")] = [realfdR; realfdR[end]]
end
realdS4df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

tmpdf = DataFrame(θ=abs.(θs))
for i = 1:4
    tmpdf[!,Symbol("dS22_$i")] = getfield.(Tuple.(imagdS4s), i)
end
for i = 1:4
    imagfdR = finitediff(getfield.(Tuple.(imagS4s), i), h)
    tmpdf[!,Symbol("finite_dS22_$i")] = [imagfdR; imagfdR[end]]
end
imagdS4df = DataFrame(θ=repeat(tmpdf[!,:θ], outer=8), var=stack(tmpdf[!,2:end])[!,1],
    val=stack(tmpdf[!,2:end])[!,2])

pr = plot(realS4df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pi = plot(imagS4df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
pdr = plot(realdS4df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("real(val)"), style(line_width=2pt));
pdi = plot(imagdS4df, x=:θ, y=:val, color=:var,
    yintercept=[0], Geom.hline(color="black", size=1pt),
    Guide.xticks(ticks=xticks), Scale.x_continuous(labels=x->@sprintf("%0.3f", x)),
    Geom.line, Guide.xlabel("θ (rad)"), Guide.ylabel("imag(val)"), style(line_width=2pt));
v = gridstack([pr pi; pdr pdi]);
v |> SVG("C:\\Users\\forrest\\Desktop\\S4s.svg", 12inch, 7inch)
