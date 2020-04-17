using Test
using LinearAlgebra
using Statistics
using StaticArrays
using CSV
using DataFrames
using GRPF
using NLsolve

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

const mₑ = 9.1093837015e-31  # kg
const qₑ = -1.602176634e-19  # C

# To completely profile in Juno:
# @profiler (for i = 1:1e6; LWMS.testfcn(); end)


########
# Derivatives
########
@testset "Derivatives" begin
    centereddiff(v, Δθ) = (v[3:end] - v[1:end-2])/(2Δθ)  # aligns to [2:end-1]

    # NOTE: as a finite diff, Δθ determines the accuracy
    # Setting too small will _greatly_ increase test runtime
    Δθ = 1e-2*π/180
    eas = [EigenAngle(th) for th in complex.(range(0, π/2; step=Δθ), -4π/180)]

    freq = Frequency(24e3)
    bfield = BField(50e-6, π/2, 0)
    ground = Ground(15, 0.001)
    electrons = Constituent(qₑ, mₑ,
            z -> waitprofile(z, 75, 0.32, H),
            z -> electroncollisionfrequency(z, H))


    M = LWMS.susceptibility(70e3, freq, bfield, electrons)

    function eval_wmatrix(eas, Δθ, M)
        Ws = NTuple{4,SArray{Tuple{2,2},ComplexF64,2,4}}[]
        dWs = NTuple{4,SArray{Tuple{2,2},ComplexF64,2,4}}[]
        for ea in eas
            T = LWMS.tmatrix(ea, M)
            W = LWMS.wmatrix(ea, T)
            dW = LWMS.dwmatrixdθ(ea, M, T)
            push!(Ws, W)
            push!(dWs, dW)
        end

        for i = 1:4
            # Extract each tuple matrix of Ws
            ws = [w[i] for w in Ws]
            dws = [dw[i] for dw in dWs]

            re_ws = Tuple.(real.(ws))
            im_ws = Tuple.(imag.(ws))

            re_dws = Tuple.(real.(dws))
            im_dws = Tuple.(imag.(dws))

            for j = 1:4
                re_diff = centereddiff(getfield.(re_ws, j), Δθ) - getfield.(re_dws, j)[2:end-1]
                im_diff = centereddiff(getfield.(im_ws, j), Δθ) - getfield.(im_dws, j)[2:end-1]
                @test maximum(abs.(re_diff)) < 1e-3
                @test maximum(abs.(im_diff)) < 1e-3
            end
        end

        return nothing
    end
    eval_wmatrix(eas, Δθ, M)


    function eval_sharpboundaryreflection(eas, Δθ, M)
        Ros = SMatrix{2,2,ComplexF64,4}[]
        dRs = SMatrix{2,2,ComplexF64,4}[]
        for ea in eas
            Ro = LWMS.sharpboundaryreflection(ea, M)
            RdR = LWMS.sharpboundaryreflection(ea, M, LWMS.Derivative_dθ())

            @test Ro ≈ RdR[SVector(1,2),:]

            push!(Ros, Ro)
            push!(dRs, RdR[SVector(3,4),:])
        end

        re_ros = Tuple.(real.(Ros))
        im_ros = Tuple.(imag.(Ros))

        re_drs = Tuple.(real.(dRs))
        im_drs = Tuple.(imag.(dRs))

        for j = 1:4
            re_diff = centereddiff(getfield.(re_ros, j), Δθ) - getfield.(re_drs, j)[2:end-1]
            im_diff = centereddiff(getfield.(im_ros, j), Δθ) - getfield.(im_drs, j)[2:end-1]
            @test maximum(abs.(re_diff)) < 1e-5
            @test maximum(abs.(im_diff)) < 1e-5
        end

        return nothing
    end
    eval_sharpboundaryreflection(eas, Δθ, M)


    modeparams = LWMS.ModeParameters(bfield, freq, ground, electrons)
    # NOTE: This function takes a relatively long time because it has to do 2*length(eas) integrations
    function eval_integratedreflection(eas, Δθ, modeparams)
        Ros = SMatrix{2,2,ComplexF64,4}[]
        dRs = SMatrix{2,2,ComplexF64,4}[]
        println("Integrating dR/dz twice for each eigenangle. This may take a minute.")
        for ea in eas
            Ro = LWMS.integratedreflection(ea, modeparams)
            RdR = LWMS.integratedreflection(ea, modeparams, LWMS.Derivative_dθ())

            @test Ro ≈ RdR[SVector(1,2),:] rtol=1e-3

            push!(Ros, Ro)
            push!(dRs, RdR[SVector(3,4),:])
        end
        println("Integrations of dR/dz complete.")

        re_ros = Tuple.(real.(Ros))
        im_ros = Tuple.(imag.(Ros))

        re_drs = Tuple.(real.(dRs))
        im_drs = Tuple.(imag.(dRs))

        for j = 1:4
            # XXX: Tolerance is so low because R we might be passing through a branch point and
            # R is unreasonably large (and probably) not defined
            # TODO: Handle branch point calculation of R
            @test centereddiff(getfield.(re_ros, j), Δθ) ≈ getfield.(re_drs, j)[2:end-1] rtol=0.1
            @test centereddiff(getfield.(im_ros, j), Δθ) ≈ getfield.(im_drs, j)[2:end-1] rtol=0.1
        end

        return nothing
    end
    eval_integratedreflection(eas, Δθ, modeparams)

    # TODO: Test for this (low priority). Tested as part of integrated reflection
    RdR = LWMS.sharpboundaryreflection(eas[1], M, LWMS.Derivative_dθ())
    LWMS.dRdθdz(RdR, (eas[1], modeparams), 75e3)


    function eval_fresnelreflection(eas, Δθ, ground, freq)
        Rgos = SMatrix{2,2,ComplexF64,4}[]
        dRgs = SMatrix{2,2,ComplexF64,4}[]
        for ea in eas
            Rgo = LWMS.fresnelreflection(ea, ground, freq)
            Rg, dRg = LWMS.fresnelreflection(ea, ground, freq, LWMS.Derivative_dθ())

            @test Rgo ≈ Rg

            push!(Rgos, Rgo)
            push!(dRgs, dRg)
        end

        re_rgos = Tuple.(real.(Rgos))
        im_rgos = Tuple.(imag.(Rgos))

        re_drgs = Tuple.(real.(dRgs))
        im_drgs = Tuple.(imag.(dRgs))

        for j = 1:4
            re_diff = centereddiff(getfield.(re_rgos, j), Δθ) - getfield.(re_drgs, j)[2:end-1]
            im_diff = centereddiff(getfield.(im_rgos, j), Δθ) - getfield.(im_drgs, j)[2:end-1]
            @test maximum(abs.(re_diff)) < 1e-4
            @test maximum(abs.(im_diff)) < 1e-4
        end

        return nothing
    end
    eval_fresnelreflection(eas, Δθ, ground, freq)


    function eval_modalequation(eas, Δθ, modeparams, ground, freq)
        Fs = ComplexF64[]
        dFs = ComplexF64[]
        for ea in eas
            RdR = LWMS.integratedreflection(ea, modeparams, LWMS.Derivative_dθ())
            Rg, dRg = LWMS.fresnelreflection(ea, ground, freq, LWMS.Derivative_dθ())

            R = RdR[SVector(1,2),:]
            dR = RdR[SVector(3,4),:]

            F = LWMS.modalequation(R, Rg)
            dF = LWMS.modalequationdθ(R, dR, Rg, dRg)

            push!(Fs, F)
            push!(dFs, dF)
        end
        # XXX: Suffering from same branch point (?) issue as `eval_integratedreflection`
        @test centereddiff(real.(Fs), Δθ) ≈ real.(dFs)[2:end-1] rtol=0.05
        @test centereddiff(imag.(Fs), Δθ) ≈ imag.(dFs)[2:end-1] rtol=0.05

        return nothing
    end
    eval_modalequation(eas, Δθ, modeparams, ground, freq)

    #
    # Test dFdθ at resonant EA explicitly
    #
    tx = Transmitter{LWMS.VerticalDipole}("", 0, 0, 0, LWMS.VerticalDipole(), LWMS.Frequency(24e3), 100e3)
    ground = LWMS.Ground(15, 0.001)
    bfield = LWMS.BField(50e-6, deg2rad(90), 0)
    electrons = Constituent(qₑ, mₑ,
            z -> waitprofile(z, 75, 0.32, H),
            z -> electroncollisionfrequency(z, H))
    modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

    # known solution w/ tolerance=1e-10
    ea = LWMS.EigenAngle(complex(1.48962303345607, -0.0127254395865))
    dFdθ, R, Rg = LWMS.solvemodalequationdθ(ea, modeparams)

    dθ = 1e-6*π/180
    Fm = LWMS.solvemodalequation(LWMS.EigenAngle(ea.θ-dθ/2), modeparams)
    Fp = LWMS.solvemodalequation(LWMS.EigenAngle(ea.θ+dθ/2), modeparams)

    finitedF = (Fp - Fm)/dθ

    # What is driving the surprisingly poor inaccuracy here?
    @test dFdθ ≈ finitedF rtol=1e-2
end

########
# No B-field
########

# TODO

# Off-diagonal terms should be 0 with no B field
@test_skip +(M[1,2], M[1,3], M[2,1], M[2,3], M[3,1], M[3,2]) == 0

########
# Vertical B-field
########
# @testset "Vertical B-field" begin
#
# Setup scenario
#
ea = EigenAngle(deg2rad(complex(85.0, -1.0)))
@test isbits(ea)
tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
# `tx` isn't just bits
ground = Ground(15, 0.001)
@test isbits(ground)
bfield = BField(50e-6, deg2rad(90), 0)
@test isbits(bfield)
electrons = Constituent(qₑ, mₑ,
        z -> waitprofile(z, 75, 0.32, H),
        z -> electroncollisionfrequency(z, H))
@test isbits(electrons) # surprisingly...

#
# susceptibility
#
M = LWMS.susceptibility(80e3, tx.frequency, bfield, electrons)

#
# tmatrix
#
T = LWMS.tmatrix(ea, M)

#
# wmatrix
#
# Make sure "analytical" solution for W matches numerical
# Tdense = Array(T)
C = ea.cosθ
L = @SMatrix [C 0 -C 0;
              0 -1 0 -1;
              0 -C 0 C;
              1 0 1 0]

W = LWMS.wmatrix(ea, T)

@test [W[1] W[3]; W[2] W[4]] ≈ 2*(L\T)*L

LWMS.bookerquartic!(ea, M)

q, B = LWMS.BOOKER_QUARTIC_ROOTS, LWMS.BOOKER_QUARTIC_COEFFS

# Confirm roots of Booker quartic satisfy ``det(Γ² + I + M) = 0``
S = ea.sinθ
S², C² = ea.sin²θ, ea.cos²θ
for i in eachindex(q)
    G = [1-q[i]^2+M[1,1] M[1,2] S*q[i]+M[1,3];
         M[2,1] 1-q[i]^2-S²+M[2,2] M[2,3];
         S*q[i]+M[3,1] M[3,2] C²+M[3,3]]
    @test isapprox(det(G), 0, atol=sqrt(eps()))
end

# Confirm Booker quartic is directly satisfied
for i in eachindex(q)
    booker = B[5]*q[i]^4 + B[4]*q[i]^3 + B[3]*q[i]^2 + B[2]*q[i] + B[1]
    @test isapprox(booker, 0, atol=sqrt(eps()))
end

#
# Sharply bounded solution
#

M = LWMS.susceptibility(95e3, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)
W = LWMS.wmatrix(ea, T)
LWMS.bookerquartic!(ea, M)

tmp = LWMS._sharpboundaryreflection(ea, M)
initR = LWMS.sharpboundaryreflection(ea, M)

@test initR[1,2] ≈ initR[2,1]

function iterativesharpR!(f, R, W)
    f .= W[2] + W[4]*R - R*W[1] - R*W[3]*R
end
initR0 = complex([1 0.1; 0.1 1])
res = nlsolve((f,x)->iterativesharpR!(f,x,W), initR0)

@test initR ≈ res.zero

# Compare to numerical solution
# BUG: This isn't working? Possibly foundational math error
function numericalsharpR(ea, q)
    sort!(q, by=LWMS.upgoing)

    S, C = ea.sinθ, ea.cosθ
    S², C² = ea.sin²θ, ea.cos²θ
    esum = zeros(ComplexF64,4)
    for i in 1:2
        Γ = [0 -q[i] 0;
             q[i] 0 -S;
             0 S 0]
        G = Γ^2 + I + M
        E = nullspace(G)
        ℋ = Γ*E  # Γℋ == -(I + M)*E

        esum += [E[1]; E[2]; ℋ[1]; ℋ[2]]
    end
    A = [C 0 -C 0; 0 1 0 1; 0 -C 0 C; 1 0 1 0]
    e = A\esum

    R = [e[3]/e[1] e[4]/e[1];
         e[3]/e[2] e[4]/e[2]]

    return R
end
@test_skip initR ≈ numericalsharpR(ea, q)

#
# Integrated ionosphere R
#
modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)
@test isbits(modeparams)

dR = LWMS.dRdz(initR, (ea, modeparams), 95e3)

R = LWMS.integratedreflection(ea, modeparams)

@test R[1,2] ≈ R[2,1]

#
# Ground reflection
#
tmp = LWMS._fresnelreflection(ea, ground, tx.frequency)
Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)

pec_ground = LWMS.Ground(1, 1e12)
vertical_ea = LWMS.EigenAngle(0)  # not necessary for PEC test?

@test abs.(LWMS.fresnelreflection(vertical_ea, pec_ground, tx.frequency)) ≈ I

#
# Modal equation
#
f = LWMS.solvemodalequation(ea, modeparams)

# known solution w/ tolerance=1e-10
# ea = LWMS.EigenAngle(complex(1.48962303345607, -0.0127254395865))
# @test isapprox(LWMS.solvemodalequation(ea, modeparams), complex(0), atol=1e-3)

#
# Find modes
#

# TODO: Make a `triangulardomain` for this problem o avoid the low right
origcoords = rectangulardomain(complex(30, -10.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)
tolerance = 1e-6

# XXX: GRPF overhead is small (<10%?) compared to time spent integrating for `R`
modes = LWMS.findmodes(origcoords, modeparams, tolerance)

for m in modes
    @test isapprox(LWMS.solvemodalequation(m, modeparams), complex(0), rtol=1e-3, atol=1e-3)
end

#
#
#

ea = modes[1]
efc = LWMS.excitationfactorconstants(ea, R, Rg, tx.frequency, ground)
hgc = LWMS.heightgains(70e3, ea, tx.frequency, efc)

dFdθ, R, Rg = LWMS.solvemodalequationdθ(ea, modeparams)
ef = LWMS.excitationfactor(ea, dFdθ, R, Rg, efc, LWMS.FC_Ez)

rx = LWMS.GroundSampler(10e3:10e3:5000e3, LWMS.FC_Ez)
E, a, p = LWMS.Efield(1000e3, modes, modeparams, tx, rx)

# XXX: Profile - this is still ridiculously slow
Ecom, phase, amp = LWMS.Efield(modes, modeparams, tx, rx)
# end






tx = Transmitter{VerticalDipole}("", 0, 0, 0, VerticalDipole(), Frequency(24e3), 100e3)
ground = Ground(15, 0.001)
bfield = BField(50e-6, deg2rad(90), 0)
electrons = Constituent(qₑ, mₑ,
        z -> waitprofile(z, 75, 0.32, H),
        z -> electroncollisionfrequency(z, H))

modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

# TODO: Make a `triangulardomain` for this problem o avoid the low right
# origcoords = rectangulardomain(complex(30, -10.0), complex(89.9, 0.0), 0.5)
origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)
tolerance = 1e-6

# XXX: GRPF overhead is small (<10%?) compared to time spent integrating for `R`
modes = LWMS.findmodes(origcoords, modeparams, tolerance)


using GRUtils



basepath = "/home/forrest/research/LAIR/ModeSolver"
# basepath = "C:\\Users\\forrest\\Desktop"
# basepath = "C:\\users\\forrest\\research\\LAIR\\ModeSolver"

raw = CSV.File(joinpath(basepath, "verticalb_day.log");
               skipto=65, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

rx = GroundSampler(dat.dist*1000, LWMS.FC_Ez)
Ecom, phase, amp = LWMS.Efield(modes, modeparams, tx, rx)

widedf = DataFrame(dist=dat.dist,
    lwpc_amp=dat.amp, lwpc_phase=dat.phase,
    lwms_amp=amp, lwms_phase=rad2deg.(phase))


plot(dat.dist, rad2deg.(phase))
hold(true)
plot(dat.dist, dat.phase)

Figure()
plot(dat.dist, dat.phase-rad2deg.(phase), yticks=(1,2))
ylim(-5, 5)



CSV.write(joinpath(basepath, "vertical.csv"), widedf)

outpath = "C:\\Users\\forrest\\UCB\\SP_2020\\PropagationModeling\\figures"

function plot(bfield, ground, tx)
    f = tx.frequency.f/1000
    B = Int(bfield.B/1e-9)
    dip = Int(rad2deg(LWMS.dip(bfield)))
    az = Int(rad2deg(LWMS.azimuth(bfield)))
    epsr = ground.ϵᵣ
    sigma = ground.σ

    gp_title = """
    TITLE = '"$f kHz\\n\\
    |B|: $B nT, dip: $(dip)°, az: $(az)°\\n\\
    h\'\': 75 km, β: 0.32\\n\\
    ϵ_r: $epsr, σ: $sigma"'"""

    open("gp_title", "w") do io
        write(io, gp_title)
    end

    ga = `gnuplot -c "$(joinpath(outpath,"vertical_amp.gp"))" "$(joinpath(basepath,"vertical.csv"))"`
    run(ga)

    gp = `gnuplot -c "$(joinpath(outpath,"vertical_phase.gp"))" "$(joinpath(basepath,"vertical.csv"))"`
    run(gp)
end
plot(bfield, ground, tx)
