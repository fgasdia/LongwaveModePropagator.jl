using LongwaveModePropagator
using LongwaveModePropagator: QE, ME
const LMP = LongwaveModePropagator

using Plots

function eyhy()
    ea=EigenAngle(1.5 - 0.1im)
    bfield=BField(50e-6, π/2, 0)
    species=Species(QE, ME,
                    z->waitprofile(z, 75, 0.32; cutoff_low=40e3),
                    electroncollisionfrequency)
    ground=Ground(15, 0.001)
    tx=Transmitter(24e3)
    rx=GroundSampler(0:5e3:2000e3, Fields.Ez)

    wvg = HomogeneousWaveguide(bfield, species, ground)
    modeequation = PhysicalModeEquation(ea, tx.frequency, wvg)
    modes = findmodes(modeequation, LMP.defaultmesh(tx.frequency))
    
    
    m₀ = LMP.referencetoground(modes[1])
    me = PhysicalModeEquation(m₀, tx.frequency, wvg)
    dFdθ, R, Rg = LMP.solvedmodalequation(me)

    t1 = abs2(1 - R[1,1]*Rg[1,1])
    t2 = abs2(1 - R[2,2]*Rg[2,2])
    t1 < t2
    eyhy1 = (1 + Rg[2,2])*R[1,2]*Rg[1,1]/((1 + Rg[1,1])*(1 - R[2,2]*Rg[2,2]))
    eyhy2 = (1 + Rg[2,2])*(1 - R[1,1]*Rg[1,1])/((1 + Rg[1,1])*R[2,1]*Rg[2,2])
    eyhy1 ≈ eyhy2

    ef = LMP.excitationfactorconstants(m₀, R, Rg, tx.frequency, ground)

    ## Heights gains
    earthradius = LMPParams().earthradius
    α = 2/earthradius
    
    z = 0
    C² = m₀.cos²θ

    expz = exp(z/earthradius)  # assumes reflection coefficients are referenced to `d = 0`

    qz = LMP.pow23(tx.frequency.k/α)*(C² + α*z)  # (k/α)^(2/3)*(C² + α*z)

    h₁z, h₂z, dh₁z, dh₂z = LMP.modifiedhankel(qz)

    F₁h₁z = ef.F₁*h₁z
    F₂h₂z = ef.F₂*h₂z

    fz = expz*(F₁h₁z + F₂h₂z)
    fy = (ef.F₃*h₁z + ef.F₄*h₂z)
    # f₂ = 1/(1im*k) df₁/dz
    fx = expz/(1im*tx.frequency.k*earthradius)*(F₁h₁z + F₂h₂z + earthradius*(ef.F₁*dh₁z + ef.F₂*dh₂z))
end
