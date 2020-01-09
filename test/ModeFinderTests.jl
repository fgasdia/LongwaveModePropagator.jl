using Test
using StaticArrays
using GRPF

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

using Gadfly
using CSV
using DataFrames

set_default_plot_size(8inch, 6inch)
font = "Segoe UI Symbol"
basetheme = Theme(major_label_font=font, major_label_font_size=14pt,
                  minor_label_font=font, minor_label_font_size=12pt,
                  key_label_font=font, key_label_font_size=12pt,
                  line_width=2pt)
Gadfly.push_theme(basetheme)

basepath = "C:\\Users\\forrest\\Desktop"

@testset "Integration Through Ionosphere" begin
    θ = deg2rad(complex(78.2520447, -3.5052794))
    freq = 24e3

    ea = LWMS.EigenAngle(θ)
    ground = LWMS.Ground(15, 0.003)
    tx = LWMS.Source(freq)

    Bmag = 0.5172267e-4
    dcl = -0.2664399
    dcm = -0.2850476
    dcn = -0.9207376
    bfield = LWMS.BField(Bmag, dcl, dcm, dcn)

    bfield = LWMS.BField(50_000e-9, 0.0, 0.0, -1.0)

    mₑ = 9.1093837015e-31  # kg
    qₑ = -1.602176634e-19  # C
    electrons = Constituent(qₑ, mₑ,
                            h -> waitprofile(h, 75, 0.32), electroncollisionfrequency)

    topheight = 92e3
    bottomheight = 0.0

    sol = LWMS.integratedreflection(ea, topheight, bottomheight, bfield, electrons, tx)

    #==
    Breaking out integratedreflection into pieces
    ==#
    Mtop = susceptibility(topheight, tx.ω, bfield, electrons)

    params = (ω=tx.ω, k=tx.k, ea=ea, species=electrons, bfield=bfield)

    # called by `sharplyboundedreflection` and dominates its runtime
    jnk = _common_sharplyboundedreflection(ea, Mtop)

    # responsible for ~1/2 of _common_sharplyboundedreflection, probably because of creation of
    # MArray and call to `roots!`
    q, B = bookerquartic(ea, Mtop)

    # called by dRdz
    T = tmatrix(ea, Mtop)
    S = smatrix(ea, T)

    Rtop = sharplyboundedreflection(ea, Mtop)

    jnk = dRdz(Rtop, params, topheight)

    prob = ODEProblem{false}(dRdz, Rtop, (topheight, bottomheight), params)
    sol = solve(prob, Tsit5(), save_everystep=false, save_start=false)
end


function unwrap!(X)
    @inbounds for i in 2:length(X)
        d = X[i] - X[i-1]
        d = d > π ? d - 2π : (d < -π ? d + 2π : d)
        X[i] = X[i-1] + d
    end
    nothing
end

#==
Vertical B field
==#
bfield = LWMS.BField(50_000e-9, 0.0, 0.0, -1.0)
tx = LWMS.Source(24e3)
rx = LWMS.Receiver("",0.0,0.0,0.0)
ground = LWMS.Ground(15, 0.001)

mₑ = 9.1093837015e-31  # kg
qₑ = -1.602176634e-19  # C
electrons = Constituent(qₑ, mₑ,
                        h -> waitprofile(h, 75, 0.32), electroncollisionfrequency)

origcoords = rectangulardomain(complex(20., -20.), complex(90.0, 0.0), 1)
origcoords .= deg2rad.(origcoords)

integrationparams = LWMS.IntegrationParameters(0.0, 91e3, bfield, tx, ground, electrons)

modes = LWMS.findmodes(origcoords,integrationparams, 1e-6)

X = range(10e3, 5000e3; step=10e3)
E, phase, amp = LWMS.Efield(X, modes, integrationparams, rx, LWMS.Ez())

raw = CSV.File("C:\\users\\forrest\\research\\LAIR\\ModeSolver\\verticalb_day.log";
               skipto=65, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

df = DataFrame(dist=vcat(X./1000, dat.dist),
               amp=vcat(amp.+150, dat.amp),
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

p = plot(df, x=df[!,:x]/1000, y=:phase, Geom.path, Scale.color_discrete,
        Guide.ylabel("phase (deg)"), Guide.xlabel("distance (km)"),
        Scale.x_continuous(format=:plain),
        # Guide.xticks(ticks=0:30:90), #Guide.yticks(ticks=0:200:1000),
        Guide.title("24 kHz\n|B|=50e3 nT, dip=90°, az=0°\nh′: 75, β: 0.32\nϵᵣ: 15, σ: 0.003"));

p |> SVGJS(joinpath(basepath, "verticalb_phase.svg"))



#==
Dipped B field
==#
bfield = LWMS.BField(50_000e-9, deg2rad(70), deg2rad(-20))
tx = LWMS.Source(24e3)
ground = LWMS.Ground(15, 0.001)

mₑ = 9.1093837015e-31  # kg
qₑ = -1.602176634e-19  # C
electrons = Constituent(qₑ, mₑ,
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
