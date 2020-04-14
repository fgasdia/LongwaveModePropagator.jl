using Test
using LinearAlgebra
using StaticArrays

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
bfield = BField(50e-6, deg2rad(90), 0)
electrons = Constituent(qₑ, mₑ,
        z -> waitprofile(z, 75, 0.32, H),
        z -> electroncollisionfrequency(z, H))

modeparams = LWMS.ModeParameters(bfield, tx.frequency, ground, electrons)

origcoords = rectangulardomain(complex(40, -10.0), complex(89.9, 0.0), 0.5)
origcoords .= deg2rad.(origcoords)
tolerance = 1e-6

modes = LWMS.findmodes(origcoords, modeparams, tolerance)

########
# Compare Booker quartic computed with M and T
# NOTE: Runtime is dominated by `roots!` so currently no meaningful difference
# between the two

ea = modes[argmax(real(getfield.(modes, :θ)))]
z = 120e3  # "TOP HT" TODO: figure out how to confirm this ht is high enough
M = LWMS.susceptibility(z, tx.frequency, bfield, electrons)
T = LWMS.tmatrix(ea, M)

qM, BM = LWMS.bookerquartic(ea, M)

qT, BT = LWMS.bookerquartic(T)

@test sort(qM, by=real) ≈ sort(qT, by=real)
@test sort(eigvals(Array(T)), by=real) ≈ sort(qT, by=real)  # eigvals is ~7 times slower than bookerquartic

# Confirm Booker quartic is directly satisfied
for i in eachindex(qT)
    booker = qT[i]^4 + BT[4]*qT[i]^3 + BT[3]*qT[i]^2 + BT[2]*qT[i] + BT[1]
    @test isapprox(booker, 0, atol=sqrt(eps()))
end

# We choose the 2 roots corresponding to upward travelling waves as being those that lie
# close to the positive real and negative imaginary axis (315° on the complex plane)
q = qT
sort!(q, by=LWMS.upgoing)
