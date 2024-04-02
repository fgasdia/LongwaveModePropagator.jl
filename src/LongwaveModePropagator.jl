module LongwaveModePropagator

import Base: @_inline_meta, @propagate_inbounds, @_propagate_inbounds_meta
import Base: ==
using Dates, LinearAlgebra

using StaticArrays
using StaticArrays: promote_tuple_eltype, convert_ntuple
using OrdinaryDiffEq, DiffEqCallbacks
using Parameters, ProgressLogging
using JSON3, StructTypes
using Interpolations
import FunctionWrappers: FunctionWrapper

using PolynomialRoots: roots!

using RootsAndPoles, Romberg, ModifiedHankelFunctionsOfOrderOneThird
export GRPFParams

# LongwaveModePropagator.jl
export propagate

# EigenAngle.jl
export EigenAngle
export attenuation, phasevelocity, referencetoground, setea

# Geophysics.jl
export BField, Species, Fields, Ground, GROUND
export waitprofile, electroncollisionfrequency, ioncollisionfrequency

# IO.jl
export ExponentialInput, TableInput, BatchInput, BasicOutput, BatchOutput

# Samplers.jl
export Receiver, Sampler, GroundSampler

# Emitters.jl
export Transmitter, Dipole, VerticalDipole, HorizontalDipole, Frequency
export inclination, azimuth

# modefinder.jl
export findmodes, PhysicalModeEquation

# Waveguides.jl
export HomogeneousWaveguide, SegmentedWaveguide


"""
    IntegrationParams{T}

Parameters passed to `OrdinaryDiffEq.jl` during the integration of the ionosphere reflection
coefficient matrix in `modefinder.jl`.

# Fields

- `tolerance::Float64 = 1e-5`: integration `atol` and `rtol`.
- `solver::T = Vern7()`: a `DifferentialEquations.jl` solver.
- `dt::Float64 = 1.0`: height step in meters (many methods use a variable step size).
- `force_dtmin::Bool = false`: if true, continue integration when solver reaches minimum
    step size.
- `maxiters::Int = 100_000`: maximum number of iterations before stopping.
"""
@with_kw struct IntegrationParams{T}
    tolerance::Float64 = 1e-5
    solver::T = Vern7()
    dt::Float64 = 1.0
    force_dtmin::Bool = false
    maxiters::Int = 100_000
end
export IntegrationParams

"""
    LMPParams{T,T2,H <: AbstractRange{Float64}}

Parameters for the `LongwaveModePropagator` module with defaults:

- `topheight::Float64 = 110e3`: starting height for integration of the ionosphere reflection
    coefficient.
- `earthradius::Float64 = 6369e3`: Earth radius in meters.
- `earthcurvature::Bool = true`: toggle inclusion of Earth curvature in calculations. This
    is not supported by all functions.
- `curvatureheight::Float64 = 50e3`: reference height for Earth curvature in meters. At this
    height, the index of refraction is 1, and is therefore the reference height for
    eigenangles.
- `approxsusceptibility::Bool = false`: use a cubic interpolating spline representation of
    [`susceptibility`](@ref) during the integration of [`dRdz`](@ref).
- `susceptibilitysplinestep::Float64 = 10.0`: altitude step in meters used to build the
    spline representation of [`susceptibility`](@ref) if `approxsusceptibility == true`.
- `grpfparams::GRPFParams = GRPFParams(100000, 1e-5, true)`: parameters for the `GRPF`
    complex root-finding algorithm.
- `integrationparams::IntegrationParams{T} =
    IntegrationParams(solver=Vern7(), tolerance=1e-5)`:
    parameters passed to `DifferentialEquations.jl` for integration of the ionosphere
    reflection coefficient.
- `wavefieldheights::H = range(topheight, 0, length=2049)`: heights in meters at which
    wavefields will be integrated.
- `wavefieldintegrationparams::IntegrationParams{T2} =
    IntegrationParams(solver=Tsit5(), tolerance=1e-6)`:
    parameters passed to `DifferentialEquations.jl` for integration of the wavefields
    used in mode conversion. The solver cannot be lazy.
- `radiationresistancecorrection::Bool = false`: Perform a radiation resistance correction
    to the calculated field amplitude for transmitter antennas with altitude > 0.

The struct is created using `Parameters.jl` `@with_kw` and supports that package's
instantiation capabilities, e.g.:

```julia
p = LMPParams()
p2 = LMPParams(earth_radius=6370e3)
p3 = LMPParams(p2; grpf_params=GRPFParams(100000, 1e-6, true))
```

See also: [`IntegrationParams`](@ref)
"""
@with_kw struct LMPParams{T,T2,H<:AbstractRange{Float64}}
    topheight::Float64 = 110e3
    earthradius::Float64 = 6369e3  # m
    earthcurvature::Bool = true
    curvatureheight::Float64 = 50e3  # m
    approxsusceptibility::Bool = false
    susceptibilitysplinestep::Float64 = 10.0
    grpfparams::GRPFParams = GRPFParams(100000, 1e-5, true)
    integrationparams::IntegrationParams{T} = IntegrationParams()
    wavefieldheights::H = range(topheight, 0; length=2049)
    wavefieldintegrationparams::IntegrationParams{T2} =
        IntegrationParams(solver=Tsit5(), tolerance=1e-6)
    radiationresistancecorrection::Bool = false
end
export LMPParams

"""
    ModeEquation

Functions can dispatch on subtypes of `ModeEquation`, although currently only
`PhysicalModeEquation` is supported.

Future work might include the `ModifiedModeEquation` of [Morfitt1976].

# References

[Morfitt1976]: D. G. Morfitt and C. H. Shellman, “‘MODESRCH’, an improved computer program
    for obtaining ELF/VLF/LF mode constants in an Earth-ionosphere waveguide,” Naval
    Electronics Laboratory Center, San Diego, CA, NELC/IR-77T, Oct. 1976. [Online].
    Available: http://www.dtic.mil/docs/citations/ADA032573.
"""
abstract type ModeEquation end  # concrete subtypes defined in modefinder files

"Indicates derivative with respect to eigenangle ``θ``."
struct Dθ end

"""
Bottom of reflection coefficient integration. `LongwaveModePropagator.jl` assumes
`BOTTOMHEIGHT = 0.0`.
"""
const BOTTOMHEIGHT = 0.0  # WARNING: if this isn't zero, many assumptions break


#
include("utils.jl")
include("Antennas.jl")
include("Emitters.jl")  # must be before EigenAngles.jl
include("EigenAngles.jl")
include("Geophysics.jl")
include("Samplers.jl")
include("TMatrix.jl")
include("Waveguides.jl")
include("Wavefields.jl")

include("bookerquartic.jl")
include("magnetoionic.jl")
include("modeconversion.jl")
include("modefinder.jl")
include("modesum.jl")

include("IO.jl")


"""
    propagate(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler;
              modes::Union{Nothing,Vector{EigenAngle}}=nothing, mesh=nothing,
              params=LMPParams())

Compute electric field `E`, `amplitude`, and `phase` at `rx`.

Precomputed waveguide `modes` can optionally be provided as a `Vector{EigenAngle}`. By
default modes are found with [`findmodes`](@ref).

If `mesh = nothing`, use [`defaultmesh`](@ref) to generate `mesh` for the
mode finding algorithm. This is ignored if `modes` is not `nothing`.
"""
function propagate(waveguide::HomogeneousWaveguide, tx::Emitter, rx::AbstractSampler;
    modes::Union{Nothing,Vector{EigenAngle}}=nothing, mesh=nothing, unwrap=true,
    params=LMPParams())

    if isnothing(modes)
        if isnothing(mesh)
            mesh = defaultmesh(tx.frequency)
        else
            if minimum(imag(mesh)) < deg2rad(-31)
                @warn "imaginary component less than -0.5410 rad (-31°) may cause wave"*
                    "fields calculated with modified Hankel functions to overflow."
            end
        end

        modeequation = PhysicalModeEquation(tx.frequency, waveguide)
        modes = findmodes(modeequation, mesh; params=params)
    end

    EE = fieldsum(modes, waveguide, tx, rx; params=params)
   
    fc = rx.fieldcomponent
    E = permutedims(EE)[:,index(fc)]  # select field components from full matrix EE

    # Type instable, but gives the expected result
    if length(E) == 1
        E = only(E)
    end

    amplitude, phase = amplitudephase(E)

    if unwrap
        unwrap!(phase)
    end

    return E, amplitude, phase
end

"""
    propagate(waveguide::SegmentedWaveguide, tx::Emitter, rx::AbstractSampler;
              mesh=nothing, params=LMPParams())

Compute electric field `E`, `amplitude`, and `phase` at `rx` through a `SegmentedWaveguide`.

If `mesh = nothing`, use [`defaultmesh`](@ref) to generate `mesh` for the
mode finding algorithm.
"""
function propagate(waveguide::SegmentedWaveguide, tx::Emitter, rx::AbstractSampler;
    mesh=nothing, unwrap=true, params=LMPParams())

    if isnothing(mesh)
        mesh = defaultmesh(tx.frequency)
    else
        if minimum(imag(mesh)) < deg2rad(-31)
            @warn "imaginary component less than -0.5410 rad (-31°) may cause wave fields"*
                "calculated with modified Hankel functions to overflow."
        end
    end

    J = length(waveguide)  # number of waveguide segments

    heighttype = typeof(params.wavefieldheights)
    wavefields_vec = Vector{Wavefields{heighttype}}(undef, J)
    adjwavefields_vec = Vector{Wavefields{heighttype}}(undef, J)

    # Calculate wavefields and adjoint wavefields for each segment of waveguide
    for j in 1:J
        wvg = waveguide[j]

        modeequation = PhysicalModeEquation(tx.frequency, wvg)

        modes = findmodes(modeequation, mesh; params=params)

        # adjoint wavefields are wavefields through adjoint waveguide, but for same modes
        # as wavefield
        adjwvg = adjoint(wvg)

        wavefields = Wavefields(params.wavefieldheights, modes)
        adjwavefields = Wavefields(params.wavefieldheights, modes)

        calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg;
                              params=params)

        wavefields_vec[j] = wavefields
        adjwavefields_vec[j] = adjwavefields
    end

    EE = fieldsum(waveguide, wavefields_vec, adjwavefields_vec, tx, rx; params=params)

    fc = rx.fieldcomponent
    E = permutedims(EE)[:,index(fc)]  # select field components from full matrix EE

    # Type instable, but gives the expected result
    if length(E) == 1
        E = only(E)
    end

    amplitude, phase = amplitudephase(E)

    if unwrap
        unwrap!(phase)
    end

    return E, amplitude, phase
end

"""
    propagate(file::AbstractString, outfile=missing; incrementalwrite=false, append=false,
              mesh=nothing)

Run the model scenario described by `file` and save the results as `outfile`.

If `outfile = missing`, the output file name will be `\$(file)_output.json`.
"""
function propagate(file::AbstractString, outfile=missing; incrementalwrite=false,
    append=false, mesh=nothing, unwrap=true, params=LMPParams())

    ispath(file) || error("$file is not a valid file name")

    if ismissing(outfile)
        basepath = dirname(file)
        filename, _ = splitext(basename(file))

        outfile = joinpath(basepath, filename)*"_output.json"
    end

    s = parse(file)
    if incrementalwrite
        s isa BatchInput || throw(ArgumentError("incrementalwrite only supported for"*
                                                "BatchInput files"))
        output = buildrunsave(outfile, s; append=append, mesh=mesh, unwrap=unwrap, params=params)
    else
        output = buildrun(s; mesh=mesh, unwrap=unwrap, params=params)

        json_str = JSON3.write(output)

        open(outfile, "w") do f
            write(f, json_str)
        end
    end

    return output
end

end
