function homogeneouswaveguidebits()
    bfield = BField(50000e-9, 0.0, 0.0, -1.0)
    ground = Ground(15, 0.001)

    electrons = Species(QE, ME,
                        h -> waitprofile(h, 75, 0.32), electroncollisionfrequency)

    waveguide = LWMS.HomogeneousWaveguide(bfield, electrons, ground)
    return isbits(waveguide)
end

@testset "Waveguides.jl" begin
    @info "Testing Waveguides"

    @test homogeneouswaveguidebits()
end
