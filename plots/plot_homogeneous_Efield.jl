using StaticArrays
using RootsAndPoles

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

using CSV
using DataFrames

basepath = "C:\\Users\\forrest\\Desktop"

@testset "Integration Through Ionosphere" begin
    θ = deg2rad(complex(78.2520447, -3.5052794))
    freq = 24e3

    ea = LWMS.EigenAngle(θ)
    ground = LWMS.Ground(15, 0.003)
    tx = LWMS.Transmitter(freq)

    Bmag = 0.5172267e-4
    dcl = -0.2664399
    dcm = -0.2850476
    dcn = -0.9207376
    bfield = LWMS.BField(Bmag, dcl, dcm, dcn)

    # bfield = LWMS.BField(50_000e-9, 0.0, 0.0, -1.0)

    mₑ = 9.1093837015e-31  # kg
    qₑ = -1.602176634e-19  # C
    electrons = Species(qₑ, mₑ,
                            h -> waitprofile(h, 75, 0.32), electroncollisionfrequency)

    #==
    Breaking out integratedreflection into pieces
    ==#
    Mtop = LWMS.susceptibility(LWMS.TOPHEIGHT, tx.frequency, bfield, electrons)

    # called by `sharplyboundedreflection` and dominates its runtime
    jnk = LWMS._sharplyboundedreflection(ea, Mtop)

    # responsible for ~1/2 of _common_sharplyboundedreflection, probably because of creation of
    # MArray and call to `roots!`
    q, B = LWMS.bookerquartic(ea, Mtop)

    Rtop = LWMS.sharplyboundedreflection(ea, Mtop)

    # called by dRdz
    T = LWMS.tmatrix(ea, Mtop)
    S = LWMS.smatrix(ea, T)

    integrationparams = LWMS.IntegrationParameters(bfield, tx.frequency, ground, electrons)
    jnk = LWMS.dRdz(Rtop, (ea, integrationparams), LWMS.TOPHEIGHT)

    prob = LWMS.ODEProblem{false}(LWMS.dRdz, Rtop, (LWMS.TOPHEIGHT, LWMS.BOTTOMHEIGHT), (ea, integrationparams))
    sol = LWMS.solve(prob, LWMS.Tsit5(), reltol=1e-8, abstol=1e-8, save_start=false)

    sol = LWMS.integratedreflection(ea, integrationparams)

    Rg = LWMS.fresnelreflection(ea, ground, tx.frequency)
    jnk = LWMS.solvemodalequation(ea, integrationparams)
end

#==
Vertical B field
==#
bfield = BField(50_000e-9, 0.0, 0.0, -1.0)
# tx = LWMS.Transmitter(24e3)
tx = Transmitter("NAA", 44.646, -67.281, 0.0, VerticalDipole(), Frequency(24e3), 100e3)
rx = GroundSampler(10e3:10e3:5000e3, LWMS.FC_Ez)
ground = Ground(15, 0.001)

θ = deg2rad(complex(78.2520447, -3.5052794))
ea = LWMS.EigenAngle(θ)

mₑ = 9.1093837015e-31  # kg
qₑ = -1.602176634e-19  # C

electrons = Species(qₑ, mₑ,
                    h -> waitprofile(h, 75, 0.32, 50e3), electroncollisionfrequency)

# Based on lwp_input.for
origcoords = rectangulardomain(complex(20, -20.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)

# angletype = eltype(origcoords)
# zroots, zpoles = grpf(θ->solvemodalequation(EigenAngle{angletype}(θ), integrationparams),
#                   origcoords, tolerance, 15000)

waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

modes = LWMS.findmodes(origcoords, tx.frequency, waveguide)

modes = [LWMS.EigenAngle(th) for th in modes]

@testset "Electric field" begin
    talt = altitude(tx)
    ralt = altitude(rx)
    rxfc = fieldcomponent(rx)

    # Transmit dipole antenna orientation with respect to propagation direction
    # See Morfitt 1980 pg 22
    Sγ, Cγ = sincos(elevation(tx))
    Sϕ, Cϕ = sincos(azimuth(tx))

    ea = modes[1]

    S = ea.sinθ  # Morfitt 1980 pg 19 specifies S is the sine of θ at `H`

    hgc = LWMS.heightgainconstants(ea, ground, tx.frequency)

    jnk = LWMS.fresnelreflection(ea, ground, tx.frequency, LWMS.Derivative_dθ())

    dFdθ, R, Rg = LWMS.solvemodalequationdθ(ea, integrationparams)

    λ₁, λ₂, λ₃, f = LWMS.excitationfactors(ea, R, Rg, dFdθ, hgc, rxfc)
    f₁, f₂, f₃ = LWMS.heightgains(talt, ea, tx.frequency, f, hgc)  # at transmitter

    # All 3 directions needed for xmtr, only "wanted" term needed for rcvr
    xmtrterm = λ₁*f₁*Cγ + λ₂*f₂*Sγ*Cϕ + λ₃*f₃*Sγ*Sϕ
    rcvrterm = LWMS.heightgains(ralt, ea, tx.frequency, f, hgc, rxfc)  # at receiver

    # NOTE: Morfitt 1980 eq 39 is different from Pappert et al 1983 which is also
    # different from Pappert and Ferguson 1986
    modesum += xmtrterm*rcvrterm*cis(-tx.frequency.k*(S - 1)*1000e3)
    end
end

E, phase, amp = LWMS.Efield(1500e3, modes, waveguide, tx, rx)

E, phase, amp = LWMS.Efield(modes, waveguide, tx, rx)

raw = CSV.File("C:\\users\\forrest\\research\\LAIR\\ModeSolver\\verticalb_day.log";
               skipto=65, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

X = rx.distance
df = DataFrame(dist=vcat(X./1000, dat.dist),
               amp=vcat(amp, dat.amp),
               phase=vcat(rad2deg.(phase), dat.phase),
               model=vcat(fill("LWMS", length(X)), fill("LWPC", length(dat.dist))))

# colors = Scale.default_discrete_colors(2);
colors = ["deepskyblue", "orange"]
p = plot(df, x=:dist, y=:amp, color=:model, Geom.path,
        Guide.ylabel("amp (dBu)"), Guide.xlabel("distance (km)"),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(format=:plain),
        Coord.cartesian(xmax=5e3),
        # Guide.xticks(ticks=0:30:90), #Guide.yticks(ticks=0:200:1000),
        Guide.title("24 kHz\n|B|=50e3 nT, dip=90°, az=0°\nh′: 75, β: 0.32\nϵᵣ: 15, σ: 0.001"));

p |> SVGJS(joinpath(basepath, "verticalb_amp.svg"))

p = plot(df, x=:dist, y=:phase, color=:model, Geom.path, Scale.color_discrete,
        Guide.ylabel("phase (deg)"), Guide.xlabel("distance (km)"),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(format=:plain),
        Coord.cartesian(xmax=5e3),
        # Guide.xticks(ticks=0:30:90), #Guide.yticks(ticks=0:200:1000),
        Guide.title("24 kHz\n|B|=50e3 nT, dip=90°, az=0°\nh′: 75, β: 0.32\nϵᵣ: 15, σ: 0.001"));

p |> SVGJS(joinpath(basepath, "verticalb_phase.svg"))



#==
Dipped B field
==#
bfield = LWMS.BField(50_000e-9, deg2rad(70), deg2rad(-20))
tx = LWMS.Source(24e3)
ground = LWMS.Ground(15, 0.001)

mₑ = 9.1093837015e-31  # kg
qₑ = -1.602176634e-19  # C
electrons = Species(qₑ, mₑ,
                        h -> waitprofile(h, 75, 0.32), electroncollisionfrequency)

origcoords = rectangulardomain(complex(20., -20.), complex(90.0, 0.0), 1)
origcoords .= deg2rad.(origcoords)

integrationparams = LWMS.IntegrationParameters(0.0, 91e3, bfield, tx, ground, electrons)

modes = LWMS.findmodes(origcoords,integrationparams, 1e-6)

X = range(10e3, 5000e3; step=10e3)
E, phase, amp = LWMS.Efield(X, modes, integrationparams)


bfield = LWMS.BField(50_000e-9, deg2rad(70), deg2rad(20))
integrationparams = LWMS.IntegrationParameters(0.0, 91e3, bfield, tx, ground, electrons)
modes = LWMS.findmodes(origcoords,integrationparams, 1e-6)
oE, ophase, oamp = LWMS.Efield(X, modes, integrationparams)


raw = CSV.File("C:\\users\\forrest\\research\\LAIR\\ModeSolver\\dippedb_day.log";
               skipto=65, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

df = DataFrame(dist=vcat(X./1000, X./1000, dat.dist),
               amp=vcat(amp.+150, oamp.+150, dat.amp),
               phase=vcat(rad2deg.(phase), rad2deg.(ophase), dat.phase),
               model=vcat(fill("LWMS", length(X)), fill("oLWMS", length(X)), fill("LWPC", length(dat.dist))))

# colors = Scale.default_discrete_colors(2);
colors = ["deepskyblue", "orange", "purple"]
p = plot(df, x=:dist, y=:amp, color=:model, Geom.path,
        Guide.ylabel("amp (dB)"), Guide.xlabel("distance (km)"),
        Scale.color_discrete_manual(colors...),
        Scale.x_continuous(format=:plain),
        # Guide.xticks(ticks=0:30:90), #Guide.yticks(ticks=0:200:1000),
        Guide.title("24 kHz\n|B|=50e3 nT, dip=70°, az=-20°\nh′: 75, β: 0.32\nϵᵣ: 15, σ: 0.001"));

p |> SVGJS(joinpath(basepath, "dippedb_amp.svg"))

p = plot(df, x=df[!,:x]/1000, y=:phase, Geom.path, Scale.color_discrete,
        Guide.ylabel("phase (deg)"), Guide.xlabel("distance (km)"),
        Scale.x_continuous(format=:plain),
        # Guide.xticks(ticks=0:30:90), #Guide.yticks(ticks=0:200:1000),
        Guide.title("24 kHz\n|B|=50e3 nT, dip=90°, az=0°\nh′: 75, β: 0.32\nϵᵣ: 15, σ: 0.003"));

p |> SVGJS(joinpath(basepath, "verticalb_phase.svg"))


#==
LWPC search region criteria (see lwp_input.for)

No idea where this actual comes from. It was added in January 1993

if f < 20 kHz
    ranger = (60, 90)
    rangei = (0, -9)
else
    ranger = (40, 90)
    rangei = (0, -6)
end
==#
