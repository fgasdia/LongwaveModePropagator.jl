"""
Map Roe values and evaluate interpolation quality.
"""

using StaticArrays
using Interpolations
using Plots

include("../src/ModeFinder.jl")
import ModeFinder
using LWMS
using IonosphereProfile

mode = ModeFinder.Mode()
mode.ω = 20e3
mode.wavenumber = mode.ω/speedoflight*1e3

drcs = ModeFinder.DirectionCosines(-0.2664399, -0.2850476, -0.9207376, 0.0)
B = 0.5172267e-4

inputs = ModeFinder.Inputs()
inputs.topheight = 95.0
inputs.bottomheight = 45.0
referenceheight = 50.0
heightinterp = 72.5

electrons = Constituent(-fundamentalcharge, mₑ,
                        h -> waitprofile(h, 75, 0.3), collisionprofile)

ζₚ = @MVector [complex(0.0), complex(0.0)]
Γₚ = @MVector [complex(0.0), complex(0.0)]
xseq = ModeFinder.XSequence(ζₚ, Γₚ)


function evalroe!(Roegrid, θ, mode::ModeFinder.Mode, inputs::LWMS.Inputs, referenceheight, heightinterp,
                 spec::IonosphereProfile.Constituent, B, drcs::ModeFinder.DirectionCosines, xseq)
    bottomheight = inputs.bottomheight
    Eud = @MArray zeros(ComplexF64, 2,2,2)
    Xhintrp = @MMatrix zeros(ComplexF64, 2, 2)

    for i in eachindex(θ)
        mode.θ = θ[i]
        X = ModeFinder.integratethroughionosphere(mode, inputs, referenceheight, spec, B, drcs)[end]
        ModeFinder.integratethroughfreespace!(X, mode, bottomheight, heightinterp, referenceheight)
        ModeFinder.eigenmatrix!(Eud, xseq, mode, false, heightinterp, referenceheight, drcs)
        R, Roe = ModeFinder.roematrix(θ[i], X, Eud)

        Roegrid[:,:,i] = log.(Roe)
    end
end


θ = [complex(i, j) for i=-90:5:90, j=-40:5:40]
Roegrid = rand(ComplexF64, 2, 2, length(θ))

evalroe!(Roegrid, θ, mode, inputs, referenceheight, heightinterp, electrons, B, drcs, xseq)

itp11real = interpolate(reshape(real(Roegrid[1,1,:]), size(θ)),  BSpline(Quadratic(Flat(OnCell()))))
v11 = itp11real(range(1,stop=19,length=1001), range(1,stop=9,length=201))
plot(contourf(v11))

plot(contourf(reshape(real(Roegrid[1,1,:]), size(θ))))

θfine = [complex(i, j) for i=-90:1:90, j=-40:1:40]
fineRoegrid = rand(ComplexF64, 2, 2, length(θfine))

evalroe!(fineRoegrid, θfine, mode, inputs, referenceheight, heightinterp, electrons, B, drcs, xseq)

plot(contourf(reshape(imag(fineRoegrid[1,1,:]), size(θfine))))
