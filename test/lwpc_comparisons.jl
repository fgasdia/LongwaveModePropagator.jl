function test_lmp(scenario)
    @unpack tx, rx, bfield, species, ground = scenario
    waveguide = HomogeneousWaveguide(bfield, species, ground)

    # me = PhysicalModeEquation(tx.frequency, waveguide)
    # modes = findmodes(me, LMP.defaultcoordinates(tx.frequency), params=params)
    #
    # E = LMP.Efield(modes, waveguide, tx, rx, params=params)
    #
    # amplitude = Vector{Float64}(undef, length(E))
    # phase = similar(amplitude)
    # @inbounds for i in eachindex(E)
    #     e = E[i]
    #     amplitude[i] = 10log10(abs2(e))  # == 20log10(abs(E))
    #     phase[i] = angle(e)  # ranges between -π:π rad
    # end

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
        resonant_horizontal_scenario=>"resonant_horizontal.log")

    for scn in (verticalB_scenario, resonant_scenario, nonresonant_scenario,
        resonant_elevatedrx_scenario, resonant_horizontal_scenario,)

        @info "    "*files[scn]
        compare(joinpath(LWPC_PATH, files[scn]), test_lmp(scn)...)
    end
end
