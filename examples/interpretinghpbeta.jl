# # Interpreting h′ and β
# 
# This example looks at the influence of h′, β, and radio frequency on amplitude curves.
# 
# ## Background
# 
# [Wait and Spies (1964)](https://archive.org/details/characteristicso300wait) is a
# technical report which characterizes the influence of Earth-ionosphere waveguide
# parameters on the propagation of VLF radio waves (the attenuation rates, phase velocities,
# etc.).
# The report is probably most cited[^1] as the source of the exponential electron
# density profile parameterized by a height h′ (km) and slope β (km⁻¹).
# The profile is usually written as:
# ```math
# N(z) = 1.43\times 10^{13} \exp(-0.15 h') \exp((\beta - 0.15)(z - h'))
# ```
# for density ``N`` (m⁻³) and altitude ``z`` (km), but as described in
# [Morfitt (1977)](https://apps.dtic.mil/docs/citations/ADA047508) or
# [Thomson (1993)](https://doi.org/10.1016/0021-9169(93)90122-F), Wait and Spies (1964)
# does not explicitly introduce this form.
# Instead, Wait introduces a "conductivity parameter" ``\omega_r`` defined by
# ```math
# \omega_r = \omega_p^2/\nu
# ```
# where ``\omega_p`` is the angular plasma frequency of the electrons and ``\nu`` is the
# effective collision frequency.
# The height at which ``\omega_r = \omega``, the angular radio wave frequency is
# approximately the wave reflection height. This statement is equivalent to ``X = Z`` from
# magnetoionic theory ([Ratcliffe (1959)](https://books.google.com/books/about/The_magneto_ionic_theory_and_its_applica.html?id=lRQvAAAAIAAJ))
# 
# Citing laboratory and some observational experiments, Wait assumes an analytical
# collision frequency profile
# ```math
# \nu = 5 \times 10^6 \exp(-0.15(z - 70)).
# ```
# Using daytime measurements made by the method of partial reflection, he calculates
# the approximation
# ```math
# \omega_r = (2.5\times 10^5)\exp(0.3(z - 70)).
# ````
# To fit measurements made at other times and locations, he parameterizes the function
# as
# ```math
# \omega_r = (2.5 \times 10^5)\exp(\beta(z - h'))
# ```
# where "``\beta`` is a constant and ``h'`` is a reference height" (the height at which
# ``\omega_r = 2.5\times 10^5``).
# He immediately suggests ``\beta = 0.3`` at day and ``\beta = 0.5`` at night, which are
# largely considered typical values today.
# Wait also parameterizes the collision frequency profile as
# ```math
# \nu(z) = \nu_0\exp(-az).
# ```
# Wait fixed ``a=0.15`` for his studies in the report.
# 
# [Morfitt and Shellman (1976)](http://www.dtic.mil/docs/citations/ADA032573)
# cites Wait and Spies and uses
# ```math
# \nu(z) = 1.82\times 10^{11} \exp(-0.15 z)
# ```
# from this, it is simple to derive the equation attributed to
# Wait and Spies (1964). The ``(\beta - 0.15)`` term comes from the assumption
# ``N(z) = N_0\exp(bz)`` and ``\nu(z) = \nu_0\exp(-az)`` where ``\beta = a + b`` and
# ``a = 0.15``. 

using Plots
using LongwaveModePropagator
const LMP = LongwaveModePropagator

## Constant values
bfield = BField(50e-6, π/2, 0)
ground = GROUND[5]

tx = Transmitter(20e3)
rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

# ## Varying h′
# 
# In this section we set the frequency to 20 kHz, β is set to 0.4 km⁻¹ and we vary h′.

function varyhp(hprimes)
    amps = Vector{Vector{Float64}}(undef, length(hprimes))
    for i in eachindex(hprimes)
        electrons = Species(LMP.QE, LMP.ME, z->waitprofile(z, hprimes[i], 0.4), electroncollisionfrequency)
        waveguide = HomogeneousWaveguide(bfield, electrons, ground)
        E, a, p = propagate(waveguide, tx, rx)
        amps[i] = a
    end
    return amps
end

hprimes = 72:1:78
amps = varyhp(hprimes)

p = plot()
function buildplots(p, amps)
    cmap = palette(:amp, length(hprimes)+1) # +1 allows us to use a darker lightest color

    for i in eachindex(hprimes)
        plot!(p, rx.distance/1000, amps[i],
              label=hprimes[i], color=cmap[i+1])
    end
end
buildplots(p, amps)
plot!(p, size=(600,400), ylims=(22, 95),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="h′")
#md savefig("interpreting_hprimes.png"); nothing # hide
#md ![](interpreting_hprimes.png)

# As h′ increases, the field strength modal interference pattern is displaced outwardly
# away from the transmitter.

# ## Varying β
# 
# Now we'll fix h′ at 78 km and vary β.

function varybeta(betas)
    amps = Vector{Vector{Float64}}(undef, length(betas))
    for i in eachindex(betas)
        electrons = Species(LMP.QE, LMP.ME, z->waitprofile(z, 78, betas[i]), electroncollisionfrequency)
        waveguide = HomogeneousWaveguide(bfield, electrons, ground)
        E, a, p = propagate(waveguide, tx, rx)
        amps[i] = a
    end
    return amps
end

betas = [0.3, 0.4, 0.5, 0.7, 0.9, 2.0]
amps = varybeta(betas)

p = plot()
function buildplots(p, amps)
    cmap = palette(:amp, length(betas)+1)

    for i in eachindex(betas)
        plot!(p, rx.distance/1000, amps[i],
              label=betas[i], color=cmap[i+1])
    end
end
buildplots(p, amps)
plot!(p, size=(600,400), ylims=(22, 95),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="β")
#md savefig("interpreting_hprimes.png"); nothing # hide
#md ![](interpreting_hprimes.png)

# Broadly, the signal increases with increasing β.
# Morfitt (1977) also suggests that an increase in β leads to the signal levels varying
# over a wider range in the regions of strong modal interference.

# ## Varying frequency
# 
# In this section we'll look at an ionosphere described with typical daytime
# Wait profile parameters, but we'll vary the radio frequency from 5 kHz to 50 kHz.

function varyfreq(freqs)
    electrons = Species(LMP.QE, LMP.ME, z->waitprofile(z, 75, 0.35), electroncollisionfrequency)
    waveguide = HomogeneousWaveguide(bfield, electrons, ground)

    amps = Vector{Vector{Float64}}(undef, length(freqs))
    for i in eachindex(freqs)
        tx = Transmitter(freqs[i])
        E, a, p = propagate(waveguide, tx, rx)
        amps[i] = a
    end
    return amps
end

freqs = 5e3:5e3:50e3
amps = varyfreq(freqs)

p = plot()
function buildplots(p, amps)
    # cdivs = [1, 2, 2, 3, 3, 3, 3, 4, 4, 4]
    # Ndiv = 6
    # d = floor(Int, 100/Ndiv)
    pal = palette(:rainbow, 7)
    cmap = [pal[1]; range(pal[2], pal[3], length=2); range(pal[4], pal[5], length=4);
            range(pal[6], pal[7], length=3)]

    for i in eachindex(freqs)
        fkHz = trunc(Int, freqs[i]/1000)
        λkm = trunc(LMP.C0/freqs[i]/1000, digits=1)
        plot!(p, rx.distance/1000, amps[i] .+ (10*i),
              label=string(fkHz, ",  ", λkm), color=cmap[i])
    end
end
buildplots(p, amps)
plot!(p, size=(600,400),
      xlabel="Range (km)", ylabel="Amplitude (dB)", legendtitle="f kHz, λ km")
#md savefig("interpreting_frequencies.png"); nothing # hide
#md ![](interpreting_frequencies.png)

# The effect of increasing frequency is similar to increasing h′.
# This makes sense because in both cases the waveguide "height" in wavelengths is
# increassing.

# [^1]: As a government technical report and not a journal paper curated by a publishing
# company, it is difficult to track how many times it has been cited by others.