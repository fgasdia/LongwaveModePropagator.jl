
function mc_scenario()
    waveguide = LMP.SegmentedWaveguide(HomogeneousWaveguide)

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          Species(QE, ME,
                                                  z -> waitprofile(z, 75, 0.32),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001)))

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          Species(QE, ME,
                                                  z -> waitprofile(z, 70, 0.25),
                                                  electroncollisionfrequency),
                                          Ground(15, 0.001), 1000e3))


    tx = Transmitter(24e3)
    rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

    E, amp, phase = LMP.bpm(waveguide, tx, rx)

    return E, amp, phase
end

const QE = LMP.QE
const ME = LMP.ME

function mc_scenario()
    waveguide = LMP.SegmentedWaveguide(HomogeneousWaveguide)

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          Species(QE, ME,
                                                  z -> waitprofile(z, 84, 1.2),
                                                  electroncollisionfrequency),
                                          Ground(81, 4.0), 0.0))

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          Species(QE, ME,
                                                  z -> waitprofile(z, 84, 0.4),
                                                  electroncollisionfrequency),
                                          Ground(81, 4.0), 500e3))

    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          Species(QE, ME,
                                                  z -> waitprofile(z, 84, 1.2),
                                                  electroncollisionfrequency),
                                          Ground(81, 4.0), 600e3))


    tx = Transmitter(VerticalDipole(), Frequency(24e3), 100e3)
    rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

    E, amp, phase = LMP.bpm(waveguide, tx, rx)

    return E, amp, phase
end

function tabular_mc_scenario()
    waveguide = LMP.SegmentedWaveguide(HomogeneousWaveguide)

    alt = 0.0:100:110e3

    collision_frequency = LMP.electroncollisionfrequency.(alt)
    density = LMP.waitprofile.(alt, 84, 1.2)
    density_itp = LinearInterpolation(alt, density, extrapolation_bc=Line())
    collision_itp = LinearInterpolation(alt, collision_frequency, extrapolation_bc=Line())
    species = Species(QE, ME, density_itp, collision_itp)
    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          species,
                                          Ground(81, 4.0), 0.0))

    collision_frequency = LMP.electroncollisionfrequency.(alt)
    density = LMP.waitprofile.(alt, 84, 0.4)
    density_itp = LinearInterpolation(alt, density, extrapolation_bc=Line())
    collision_itp = LinearInterpolation(alt, collision_frequency, extrapolation_bc=Line())
    species = Species(QE, ME, density_itp, collision_itp)
    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          species,
                                          Ground(81, 4.0), 500e3))

    collision_frequency = LMP.electroncollisionfrequency.(alt)
    density = LMP.waitprofile.(alt, 84, 1.2)
    density_itp = LinearInterpolation(alt, density, extrapolation_bc=Line())
    collision_itp = LinearInterpolation(alt, collision_frequency, extrapolation_bc=Line())
    species = Species(QE, ME, density_itp, collision_itp)
    push!(waveguide, HomogeneousWaveguide(BField(50e-6, π/2, 0),
                                          species,
                                          Ground(81, 4.0), 600e3))

    tx = Transmitter(VerticalDipole(), Frequency(24e3), 100e3)
    rx = GroundSampler(0:5e3:2000e3, Fields.Ez)

    E, amp, phase = LMP.bpm(waveguide, tx, rx)

    return E, amp, phase
end
