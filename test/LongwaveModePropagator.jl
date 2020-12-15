function test_propagate_segmented(scenario)
    @unpack tx, rx, bfield, species, ground, distances = scenario

    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i], ground[i],
                                                         distances[i]) for i in 1:2])

    E, amp, phase = propagate(waveguide, tx, rx)
end

@testset "LongwaveModePropagator.jl" begin
    @info "Testing LongwaveModePropagator"

end
