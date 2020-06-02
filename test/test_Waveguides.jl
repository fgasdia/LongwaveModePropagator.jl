function homogeneouswaveguidebits()
    bfield = BField(50_000e-9, 0.0, 0.0, -1.0)
    ground = Ground(15, 0.001)

    mₑ = 9.1093837015e-31  # kg
    qₑ = -1.602176634e-19  # C

    electrons = Species(qₑ, mₑ,
                        h -> waitprofile(h, 75, 0.32), electroncollisionfrequency)

    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    return isbits(waveguide)
end

@testset "Waveguides.jl" begin
    @info "Testing Waveguides"

    @test homogeneouswaveguidebits()
end
