using Dates
using Plots
using LongwaveModePropagator, GeographicLib
const LMP = LongwaveModePropagator

using LMPTools, PropagationModelPrep, SubionosphericVLFInversionAlgorithms
const SIA = SubionosphericVLFInversionAlgorithms


function buildpaths()
    transmitters = [TRANSMITTER[:NLK], TRANSMITTER[:NML]]
    receivers = [
        Receiver("Whitehorse", 60.724, -135.043, 0.0, VerticalDipole()),
        Receiver("Churchill", 58.74, -94.085, 0.0, VerticalDipole()),
        Receiver("Stony Rapids", 59.253, -105.834, 0.0, VerticalDipole()),
        Receiver("Fort Smith", 60.006, -111.92, 0.0, VerticalDipole()),
        Receiver("Bella Bella", 52.1675508, -128.1545219, 0.0, VerticalDipole()),
        Receiver("Nahanni Butte", 61.0304412, -123.3926734, 0.0, VerticalDipole()),
        Receiver("Juneau", 58.32, -134.41, 0.0, VerticalDipole()),
        Receiver("Ketchikan", 55.35, -131.673, 0.0, VerticalDipole()),
        Receiver("Winnipeg", 49.8822, -97.1308, 0.0, VerticalDipole()),
        Receiver("IslandLake", 53.8626, -94.6658, 0.0, VerticalDipole()),
        Receiver("Gillam", 56.3477, -94.7093, 0.0, VerticalDipole())
    ]
    paths = [(tx, rx) for tx in transmitters for rx in receivers]

    return paths
end

struct LMPOutputs
    output_ranges::Vector{Float64}
    amplitude::Vector{Float64}
    phase::Vector{Float64}
end

# I can't replicate the problem anymore with homogeneous ionospheres. Maybe do a more
# rigorous test with many random ionospheres?
function homogeneous()
    bfield = BField(5.6e-5, 0.17, 0.129, -0.977)
    ground = GROUND[5]
    species = Species(LMP.QE, LMP.ME, z->waitprofile(z, 73.714, 0.2966), electroncollisionfrequency)
    wvg = SegmentedWaveguide([
        HomogeneousWaveguide(bfield, species, ground, 0.0),
        HomogeneousWaveguide(bfield, species, ground, 1e3)
    ])
    tx = TRANSMITTER[:NLK]
    gs = GroundSampler(0:5e3:2000e3, Fields.Ez)
    _, a, p = propagate(wvg, tx, gs)
    LMPOutputs(gs.distance, a, p)
end
hlmp = homogeneous()
plot(hlmp.output_ranges/1000, hlmp.amplitude)
plot!(xlims=(900,1100), ylims=(44,45))


# I can't replicate the problem anymore with homogeneous ionospheres. Maybe do a more
# rigorous test with many random ionospheres?
function groundchange()
    bfield = BField(5.6e-5, 0.17, 0.129, -0.977)
    species = Species(LMP.QE, LMP.ME, z->waitprofile(z, 73.714, 0.2966), electroncollisionfrequency)
    wvg = SegmentedWaveguide([
        HomogeneousWaveguide(bfield, species, GROUND[7], 0.0),
        HomogeneousWaveguide(bfield, species, GROUND[5], 1e3)
    ])
    tx = TRANSMITTER[:NLK]
    gs = GroundSampler(0:5e3:2000e3, Fields.Ez)
    _, a, p = propagate(wvg, tx, gs)
    LMPOutputs(gs.distance, a, p)
end
hlmp = groundchange()
plot(hlmp.output_ranges/1000, hlmp.amplitude)
plot!(xlims=(900,1100), ylims=(44,45))

# groundchange wavefields

zs = wavefields_vec[1].heights
p = plot()
for j = 1:3
    bfield = BField(5.6e-5, 0.17, 0.129, -0.977)
    species = Species(LMP.QE, LMP.ME, z->waitprofile(z, 73.714, 0.2966), electroncollisionfrequency)

    ground = GROUND[7]
    wvg = HomogeneousWaveguide(bfield, species, ground, 0.0)
    me = PhysicalModeEquation(EigenAngle(0), tx.frequency, wvg)
    ea = first(findmodes(me))
    
    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species; params=params)
    plot!(real(getindex.(e,3)), zs)
    # R = LMP.bookerreflection(ea, e[end])
    # Rg = LMP.fresnelreflection(ea, ground, tx.frequency)
    # b1, b2 = LMP.boundaryscalars(R, Rg, e[end], LMP.isisotropic(bfield))
    # @info "" R Rg b1 b2
end
p













function inp(lo, la, geoaz)
    dt = DateTime(2020, 3, 1, 20, 00)  # day
    bfield = igrf(geoaz, la, lo, year(dt))
    ground = GROUND[LMPTools.get_groundcode(la, lo)]
    h, b = ferguson(la, zenithangle(la, lo, dt), dt)
    # @info "" h b bfield ground
    species = Species(LMP.QE, LMP.ME, z->waitprofile(z, h, b), electroncollisionfrequency)

    return bfield, species, ground
end

function fulllmp()
    paths = buildpaths()
    N = 2

    dat = Vector{LMPOutputs}(undef, length(N))
    idx = 1
    for i in N
        tx, rx = paths[i]
        _, wpts = SIA.pathpts(tx, rx; dist=100e3)
        geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

        wvg = SegmentedWaveguide([
            HomogeneousWaveguide(inp(wpts[i].lon, wpts[i].lat, geoaz)..., wpts[i].dist) for i in eachindex(wpts)
        ])

        gs = GroundSampler(0:5e3:round(range(tx, rx)+10e3, digits=-4, RoundUp), Fields.Ez)
        _, a, p = propagate(wvg, tx, gs)
        dat[idx] = LMPOutputs(gs.distance, a, p)
        idx += 1
    end
    return dat
end


flmp = fulllmp()

p = plot(legend=false, xlims=(1300, 1500), ylims=(44, 50))
[plot!(o.output_ranges/1000, o.amplitude) for o in flmp]
p

p = plot(legend=false, xlims=(1300, 1500), ylims=(320, 340))
[plot!(o.output_ranges/1000, rad2deg.(o.phase)) for o in flmp]
p




####
paths = buildpaths()
tx, rx = paths[2] 
_, wpts = SIA.pathpts(tx, rx; dist=100e3)
geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

params = LMPParams(integrationparams=IntegrationParams(tolerance=1e-9),
    wavefieldheights=range(110e3,0,length=2049),
    wavefieldintegrationparams=IntegrationParams(solver=LMP.Vern7(lazy=false), tolerance=1e-10))

# All fields (including ground) change
ii = findfirst(isequal(1300e3), getfield.(wpts, :dist))

wvg = SegmentedWaveguide([
    HomogeneousWaveguide(inp(wpts[ii].lon, wpts[ii].lat, geoaz)..., 0.0),
    HomogeneousWaveguide(inp(wpts[ii+1].lon, wpts[ii+1].lat, geoaz)..., wpts[ii+1].dist),
    HomogeneousWaveguide(inp(wpts[ii+2].lon, wpts[ii+2].lat, geoaz)..., wpts[ii+2].dist)
])

gs = GroundSampler(0:500:2000e3, Fields.Ez)
E, a, p = propagate(wvg, tx, gs; params)

p = plot(legend=false, xlims=(1300, 1500), ylims=(44, 50))
plot!(gs.distance/1000, a)

me0 = PhysicalModeEquation(EigenAngle(0), tx.frequency, wvg[1])
m0 = getfield.(findmodes(me0), :θ)

me1 = PhysicalModeEquation(EigenAngle(0), tx.frequency, wvg[2])
m1 = getfield.(findmodes(me1), :θ)

me2 = PhysicalModeEquation(EigenAngle(0), tx.frequency, wvg[3])
m2 = getfield.(findmodes(me2), :θ)

rad2deg.(m2 .- m1)
rad2deg.(m1 .- m0)

# The modes are all less than 1 degree apart and there is the same number of them.
# Maybe the problem is in wavefield integration or conversion coeffs?

####
waveguide = wvg;
params = LMPParams(wavefieldintegrationparams=IntegrationParams(dt=1e-9), wavefieldheights=range(110e3,0.0,length=1025))

# From LongwaveModePropagator `propagate`
J = length(waveguide)  # number of waveguide segments

heighttype = typeof(params.wavefieldheights)
wavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, J)
adjwavefields_vec = Vector{LMP.Wavefields{heighttype}}(undef, J)

# Calculate wavefields and adjoint wavefields for each segment of waveguide
for j in 1:J
    wvg = waveguide[j]

    modeequation = PhysicalModeEquation(tx.frequency, wvg)

    modes = findmodes(modeequation; params=params)

    # adjoint wavefields are wavefields through adjoint waveguide, but for same modes
    # as wavefield
    adjwvg = LMP.adjoint(wvg)

    wavefields = LMP.Wavefields(params.wavefieldheights, modes)
    adjwavefields = LMP.Wavefields(params.wavefieldheights, modes)

    LMP.calculate_wavefields!(wavefields, adjwavefields, tx.frequency, wvg, adjwvg;
                          params=params)

    wavefields_vec[j] = wavefields
    adjwavefields_vec[j] = adjwavefields
end

j = 1
p = plot()
for i = 1:6
    plot!(real(getindex.(wavefields_vec[j].v, i)[:,1]), wavefields_vec[j].heights)
end
p
# Interestingly the wavefields appear very similar for wvgs 1/2 and different for 0/1, 0/2
####

zs = wavefields_vec[1].heights
p = plot()
for j = 1:3
    ea = wavefields_vec[j].eas[1]
    bfield = waveguide[j].bfield
    species = waveguide[j].species
    ground = waveguide[j].ground
    
    e = LMP.integratewavefields(zs, ea, tx.frequency, bfield, species; params=params)
    plot!(real(getindex.(e,3)), zs)
    # R = LMP.bookerreflection(ea, e[end])
    # Rg = LMP.fresnelreflection(ea, ground, tx.frequency)
    # b1, b2 = LMP.boundaryscalars(R, Rg, e[end], LMP.isisotropic(bfield))
    # @info "" R Rg b1 b2
end
p




#####

# From `Efield`
X = LMP.distance(gs, tx)
maxX = maximum(X)
Xlength = length(X)
E = Vector{ComplexF64}(undef, Xlength)

frequency = tx.frequency
k = frequency.k

Q = 0.6822408*sqrt(frequency.f*tx.power)

# Initialize
J = length(waveguide)
M = 0  # number of eigenangles in previous segment. Current segment is N
xmtrfields = Vector{ComplexF64}(undef, 0)  # fields generated by transmitter
previous_xmtrfields = similar(xmtrfields)  # fields saved from previous segment
rcvrfields = similar(xmtrfields)  # fields at receiver location

i = 1  # index of X
for j = 1:J  # index of waveguide
    wvg = waveguide[j]
    wavefields = wavefields_vec[j]
    eas = LMP.eigenangles(wavefields)
    N = LMP.nummodes(wavefields)

    # Identify distance at beginning of segment
    segment_start = wvg.distance
    # maxX < segment_start && break  # no farther X; break

    # Identify distance at end of segment
    if j < J
        segment_end = waveguide[j+1].distance
    else
        # last segment
        segment_end = typemax(typeof(segment_start))
    end

    # xmtrfields is for `Hy`
    resize!(xmtrfields, N)
    resize!(rcvrfields, N)
    if j > 1
        adjwavefields = adjwavefields_vec[j]
        prevwavefields = wavefields_vec[j-1]
        conversioncoeffs = LMP.modeconversion(prevwavefields, wavefields, adjwavefields;
                                          params=params)
    end

    # Calculate the mode terms (height gains and excitation factors) up to the current
    # segment
    for n = 1:N
        modeequation = PhysicalModeEquation(eas[n], frequency, wvg)
        txterm, rxterm = LMP.modeterms(modeequation, tx, gs; params=params)
        if j == 1
            # Transmitter exists only in the transmitter slab (obviously)
            xmtrfields[n] = txterm
        else
            # Otherwise, mode conversion of transmitted fields
            xmtrfields_sum = zero(eltype(xmtrfields))
            for m = 1:M
                xmtrfields_sum += previous_xmtrfields[m]*conversioncoeffs[m,n]
            end
            xmtrfields[n] = xmtrfields_sum
        end

        rcvrfields[n] = xmtrfields[n]*rxterm
    end

    # Calculate E at each distance in the current waveguide segment
    while X[i] < segment_end
        x = X[i] - segment_start
        factor = Q/sqrt(abs(sin(X[i]/params.earthradius)))

        totalfield = zero(eltype(E))
        for n = 1:N
            S₀ = referencetoground(eas[n].sinθ; params=params)
            totalfield += rcvrfields[n]*cis(-k*x*(S₀ - 1))*factor
        end

        E[i] = totalfield
        i += 1
        i > Xlength && break
    end

    # If we've reached the end of the current segment and there are more segments,
    # prepare for next segment
    if j < J
        # End of current slab
        x = segment_end - segment_start

        resize!(previous_xmtrfields, N)
        for n = 1:N
            S₀ = referencetoground(eas[n].sinθ; params=params)

            # Excitation factors at end of slab
            xmtrfields[n] *= cis(-k*x*(S₀ - 1))
            previous_xmtrfields[n] = xmtrfields[n]
        end
        M = N  # set previous number of modes
    end
end

# At transmitter (within 1 meter from it), E is complex NaN or Inf
if X[1] < 1
    # Used in LWPC `lw_sum_modes.for`, but not sure where they got it
    # amplitude = 10log10(80*Q)
    E[1] = sqrt(80*Q) + 0.0im # == 10^(amplitude/20)
end

















####


function comparemodels(lwpc)
    dt = DateTime(2020, 3, 1, 20, 00)  # day
    hbfcn(lo, la, dt) = ferguson(la, zenithangle(la, lo, dt), dt)
    
    pathstep = 100e3
    paths = buildpaths()

    batch = BatchInput{ExponentialInput}()
    batch.name = "estimate"
    batch.description = ""
    batch.datetime = Dates.now()
    batch.inputs = Vector{ExponentialInput}(undef, length(paths))

    for i in eachindex(paths)
        tx, rx = paths[i]
        input = SIA.model_observation(hbfcn, tx, rx, dt; pathstep)
        batch.inputs[i] = input
    end

    if lwpc
        computejob = LocalParallel("estimate", ".", "C:\\LWPCv21\\lwpm.exe", 16, 90)
        output = LWPC.run(batch, computejob; savefile=false)
    else
        output = LMP.buildrun(batch; params=LMPParams(approxsusceptibility=true))
    end

    return output
end

lwpco = comparemodels(true)
lmpo = comparemodels(false)

p = plot(legend=false, xlims=(1e3, 1.5e3), ylims=(40, 50), xlabel="Range", ylabel="Amplitude")
[plot!(o.output_ranges/1000, o.amplitude) for o in lmpo.outputs]

p = plot(legend=false)
[plot!(o.output_ranges, o.amplitude, color="red") for o in lwpco.outputs]

p = plot(legend=false, xlabel="Range", ylabel="Amplitude")
for i = 1:length(lmpo.outputs)
    mo = lmpo.outputs[i]
    wo = lwpco.outputs[i]
    plot!(mo.output_ranges/1000, mo.amplitude .- wo.amplitude)
end


amps = Vector{Float64}(undef, length(paths))
phases = similar(amps)
for i in eachindex(paths)
    tx, rx = paths[i]
    d = range(tx, rx)
    o = output.outputs[i]

    # Using LMP we could skip the LinearInterpolation and compute the field at exactly
    # the correct distance, but for consistency with LWPC we'll interpolate the output.
    # NOTE: step size here should match `output_ranges` step in `model_observation`!
    aitp = LinearInterpolation(0:5e3:last(o.output_ranges), o.amplitude)
    pitp = LinearInterpolation(0:5e3:last(o.output_ranges), o.phase)
    amps[i] = aitp(d)
    phases[i] = pitp(d)
end
