using StaticArrays
using RootsAndPoles

using LongwaveModeSolver
const LWMS = LongwaveModeSolver

using CSV
using DataFrames


#==
Vertical B field
==#
bfield = BField(50_000e-9, 0.0, 0.0, -1.0)
tx = Transmitter("NAA", 44.646, -67.281, 0.0, VerticalDipole(), Frequency(24e3), 100e3)
ground = Ground(15, 0.001)

mₑ = 9.1093837015e-31  # kg
qₑ = -1.602176634e-19  # C

electrons = Species(qₑ, mₑ,
                    h -> waitprofile(h, 75, 0.32, 50e3), electroncollisionfrequency)

# Based on lwp_input.for
origcoords = rectangulardomain(complex(20, -20.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)

waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)

modes = LWMS.findmodes(origcoords, tx.frequency, waveguide)


basepath = "/home/forrest/research/LAIR/ModeSolver"

raw = CSV.File(joinpath(basepath, "verticalb_day.log");
               skipto=65, delim=' ', ignorerepeated=true, header=false)

dat = DataFrame(dist=vcat(raw.Column1, raw.Column4, raw.Column7),
                amp=vcat(raw.Column2, raw.Column5, raw.Column8),
                phase=vcat(raw.Column3, raw.Column6, raw.Column9))

rx = GroundSampler(dat.dist*1000, LWMS.FC_Ez)

modes = EigenAngle.(modes)
Ecom, phase, amp = LWMS.Efield(modes, waveguide, tx, rx)

widedf = DataFrame(dist=dat.dist,
                   lwpc_amp=dat.amp, lwpc_phase=dat.phase,
                   lwms_amp=amp, lwms_phase=rad2deg.(phase))

CSV.write(joinpath(basepath, "vertical.csv"), widedf)

outpath = "/home/forrest/UCB/SP_2020/PropagationModeling/figures"

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

    ga = `gnuplot -c "$(joinpath(outpath,"vertical_amp_linux.gp"))" "$(joinpath(basepath,"vertical.csv"))" "$(joinpath(outpath,""))"`
    run(ga)

    gp = `gnuplot -c "$(joinpath(outpath,"vertical_phase_linux.gp"))" "$(joinpath(basepath,"vertical.csv"))" "$(joinpath(outpath,""))"`
    run(gp)
end
plot(bfield, ground, tx)




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

modes = LWMS.findmodes(origcoords,integrationparams, 1e-6)

X = range(10e3, 5000e3; step=10e3)
E, phase, amp = LWMS.Efield(X, modes, )


bfield = LWMS.BField(50_000e-9, deg2rad(70), deg2rad(20))
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
Ecom, phase, amp = LWMS.Efield(EigenAngle.(modes), modeparams, tx, rx)

widedf = DataFrame(dist=dat.dist,
    lwpc_amp=dat.amp, lwpc_phase=dat.phase,
    lwms_amp=amp, lwms_phase=rad2deg.(phase))


GRUtils.plot(dat.dist, rad2deg.(phase))
hold(true)
GRUtils.plot(dat.dist, dat.phase)

Figure()
GRUtils.plot(dat.dist, dat.phase-rad2deg.(phase), yticks=(1,2))
ylim(-5, 5)
