function test_lmp(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    E, amplitude, phase = propagate(waveguide, tx, rx)

    lmp_d = collect(distance(rx, tx)/1e3)
    lmp_a = amplitude
    lmp_p = rad2deg.(phase)

    return lmp_d, lmp_a, lmp_p
end

function test_lmp_segmented(scenario)
    @unpack tx, rx, bfield, species, ground, distances = scenario

    waveguide = SegmentedWaveguide([HomogeneousWaveguide(bfield[i], species[i], ground[i],
                                                         distances[i]) for i in 1:2])

    E, amplitude, phase = propagate(waveguide, tx, rx)

    lmp_d = collect(distance(rx, tx)/1e3)
    lmp_a = amplitude
    lmp_p = rad2deg.(phase)

    return lmp_d, lmp_a, lmp_p
end

function compare(lwpc_file, lmp_d, lmp_a, lmp_p)
    lwpc_d, lwpc_a, lwpc_p = readlog(lwpc_file)

    distmask = lwpc_d .> 300
    lm_p = mod.(lmp_p[distmask], 360)  # modulo because we're not interested in wrapping
    lw_p = mod.(lwpc_p[distmask], 360)

    @test lmp_d ≈ lwpc_d
    @test meanabsdiff(lmp_a[distmask], lwpc_a[distmask]) < 0.4
    @test meanabsdiff(lm_p, lw_p) < 4.0
end


@testset "LWPC comparisons" begin
    @info "Comparing to LWPC..."

    files = Dict(verticalB_scenario=>"verticalB.log", resonant_scenario=>"resonant.log",
        nonresonant_scenario=>"nonresonant.log",
        resonant_elevatedrx_scenario=>"resonant_elevatedrx.log",
        resonant_horizontal_scenario=>"resonant_horizontal.log",
        segmented_scenario=>"segmented.log")

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        resonant_elevatedrx_scenario, resonant_horizontal_scenario,)

        @info "    "*files[scn]
        compare(joinpath(LWPC_PATH, files[scn]), test_lmp(scn)...)
    end
    for scn in (segmented_scenario,)
        @info "    "*files[scn]
        compare(joinpath(LWPC_PATH, files[scn]), test_lmp_segmented(scn)...)
    end
end
