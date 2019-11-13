using Test
using StaticArrays

using LongwaveModeSolver
const LWMS = LongwaveModeSolver


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
