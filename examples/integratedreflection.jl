# # Solvers for ionosphere reflection coefficient
#
# LongwaveModePropagator.jl uses the technique presented by
# [Budden (1955)](https://doi.org/10.1098/rspa.1955.0027) to compute the reflection
# coefficient of a horiontally stratified ionosphere consisting of an
# anisotropic, collisional plasma.
#
# Budden (1955) derives a differential equation for the reflection coefficient
# ``\bm{R}``:
# ```math
# \frac{2i}{k}\frac{\mathrm{d}\bm{R}}{\mathrm{d}z} = W_{21} + W_{22}\bm{R} - \bm{R}W_{11} - \bm{R}W_{12}\bm{R}
# ```
# where ``W`` is a ``4\times 4`` matrix divided into four ``2\times 2`` submatrices
# each containing components of the ``\bm{T}`` matrix.
#
# To obtain the reflection coefficient ``\bm{R}`` at any level, we integrate
# ``\mathrm{d}\bm{R}/\mathrm{d}z`` downward through the ionosphere from a great height
# and the height at which we stop is ``\bm{R}`` for a sharp boundary at that level with
# free space below.
#
# The reflection coefficient matrix ``\bm{R}`` consists of four complex reflection
# coefficients for the different combinations of incident and reflected wave
# polarization. These reflection coefficients are also functions of the wave
# angle of incidence, and therefore the integration of ``\mathrm{d}\bm{R}/\mathrm{d}z``
# occurs tens of thousands to hundreds of thousands of times for every
# `HomogeneousWaveguide` in the GRPF mode finder. It is therefore extremely important
# that the differential equations solver be as efficient as possible to minimize
# runtime.
# The mode finder is responsible for more than 90% of the runtime of
# LongwaveModePropagator.jl, and most of the mode finder runtime is the ionosphere
# reflection coefficient integration.
#
# In this example, we will compare solvers and tolerances from
# DifferentialEquations to determine the most efficient
# combination with reasonable robustness and accuracy.
#
# ## Setup
#
# First we load the packages we need.

using Statistics
using Plots, DisplayAs

using ..LongwaveModePropagator
using ..LongwaveModePropagator: StaticArrays
using ..LongwaveModePropagator: QE, ME, integratedreflection
using ..LongwaveModePropagator: RK4, Tsit5, BS5, OwrenZen5, Vern6, Vern7, Vern8

# We will evaluate the solutions across a range of different random ionospheres,
# frequencies, and angles of incidence.
# Each scenario is described by a `PhysicalModeEquation`.

function generatescenarios(N)
    eas = EigenAngle.(complex.(rand(N)*π/2, rand(N)*deg2rad(-10)))
    frequencies = Frequency.(abs.(randn(N)*60e3))

    B = rand(30e-6:5e-7:60e-6, N)
    bfields = BField.(B, rand(N)*π/2, rand(N)*2π)

    hps = randn(N)*5 .+ 79
    betas = randn(N)*0.2 .+ 0.45

    scenarios = Vector{PhysicalModeEquation}(undef, N)
    for i = 1:N
        species = Species(QE, ME, z->waitprofile(z, hps[i], betas[i]),
                          electroncollisionfrequency)
        ground = GROUND[5]  ## not used
        waveguide = HomogeneousWaveguide(bfields[i], species, ground)

        me = PhysicalModeEquation(eas[i], frequencies[i], waveguide)
        scenarios[i] = me
    end

    return scenarios
end

scenarios = generatescenarios(30);

# ## Reference solutions
#
# To evaluate the accuracy of the reflection coefficients, we compare to a very
# low tolerance Runge-Kutta Order 4 method. The DifferentialEquations.jl implementation
# of `RK4` uses adaptive stepping.

ip = IntegrationParams(tolerance=1e-14, solver=RK4(), maxiters=1_000_000)
params = LMPParams(integrationparams=ip)

Rrefs = [integratedreflection(scenario, params=params) for scenario in scenarios];

# ## Evaluate solvers
#
# Now let's compute and time the results a set of different methods for a range
# of tolerances.
# We repeat the integration `N = 25` times for each combination of parameters
# to get a more accurate average time.

function compute(scenarios, tolerances, solvers)
    dims = length(scenarios), length(tolerances), length(solvers)
    Rs = Array{StaticArrays.SMatrix{2,2,ComplexF64,4}}(undef, dims...)
    times = Array{Float64}(undef, dims...)

    for k in eachindex(solvers)
        for j in eachindex(tolerances)
            ip = IntegrationParams(tolerance=tolerances[j], solver=solvers[k])
            params = LMPParams(integrationparams=ip)

            for i in eachindex(scenarios)
                ## warmup
                R = integratedreflection(scenarios[i], params=params)

                ## loop for average time
                N = 25
                t0 = time_ns()
                for n = 1:N
                    R = integratedreflection(scenarios[i], params=params)
                end
                ttotal = time_ns() - t0

                Rs[i,j,k] = R
                times[i,j,k] = ttotal/N
            end
        end
    end
    return Rs, times
end

tolerances = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
tolerancestrings = string.(tolerances)

solvers = [RK4(), Tsit5(), BS5(), OwrenZen5(), Vern6(), Vern7(), Vern8()]
solverstrings = replace.(string.(solvers), "OrdinaryDiffEq."=>"")

Rs, times = compute(scenarios, tolerances, solvers);

# We'll measure the error in the reflection coefficient matrices
# as the maximum absolute difference of the four elements of the matrix
# compared to the reference reflection coefficient matrix.

function differr(a, ref)
    differences = a .- ref
    aerror = similar(a, Float64)
    for i in eachindex(differences)
        absdiff = abs.(differences[i])
        aerror[i] = maximum(absdiff)
    end
    return aerror
end

Rerrs = differr(Rs, Rrefs)
mean_Rerrs = dropdims(mean(Rerrs, dims=1), dims=1)

img = heatmap(tolerancestrings, solverstrings, permutedims(log10.(mean_Rerrs)),
              clims=(-8, -3),
              xlabel="tolerance", ylabel="solver",
              colorbar_title="log₁₀ max abs difference")
DisplayAs.PNG(img)  #hide

# And the average runtimes are

mean_times = dropdims(mean(times, dims=1), dims=1)

img = heatmap(tolerancestrings, solverstrings, permutedims(mean_times)/1e6,
              clims=(0, 10),
              xlabel="tolerance", ylabel="solver", colorbar_title="time (μs)")
DisplayAs.PNG(img)  #hide

# The best combination of runtime and accuracy occurs at roughly `Vern7` with a
# tolerance of `1e-8`.
# The integration accuracy improves considerably for relatively little additional
# computation time if the tolerance is changed to `1e-9`, but that accuracy is not
# usually required.
