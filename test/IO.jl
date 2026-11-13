function generate_basic(;write=true)
    # Waveguide definition
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    hprimes = [72.0, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]
    b_mags = fill(50e-6, length(segment_ranges))
    b_dips = fill(π/2, length(segment_ranges))
    b_azs = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    frequency = 24e3

    output_ranges = collect(0:20e3:2000e3)

    input = ExponentialInput()
    input.name = "basic"
    input.description = "Test ExponentialInput"
    input.datetime = Dates.now()
    input.segment_ranges = segment_ranges
    input.hprimes = hprimes
    input.betas = betas
    input.b_mags= b_mags
    input.b_dips = b_dips
    input.b_azs = b_azs
    input.ground_sigmas = ground_sigmas
    input.ground_epsrs = ground_epsrs
    input.frequency = frequency
    input.output_ranges = output_ranges

    write && JSON.json("basic.json", input; pretty=true)

    return input
end

function generate_basic_fieldcomponent(fieldcomponent=Fields.Ez; write=true)
    input = generate_basic(write=false)
    input.fieldcomponent = fieldcomponent
    
    if fieldcomponent == Fields.Ez
        write &&JSON.json("basic.json", input; pretty=true)
    elseif fieldcomponent == Fields.E
        write && JSON.json("basic-E.json", input; pretty=true)
    elseif fieldcomponent == Fields.EH
        write && JSON.json("basic-EH.json", input; pretty=true)
    else
        error("Invalid field component: $fieldcomponent")
    end

    return input
end

function generate_basic_power(power=1000.0; write=true)
    input = generate_basic(write=false)
    input.power = power

    if power == 1000.0
        write && JSON.json("basic.json", input; pretty=true)
    else
        write && JSON.json("basic-power.json", input; pretty=true)
    end

    return input
end

function generate_table(;write=true)
    # Waveguide definition
    segment_ranges = [0, 500e3, 1000e3, 1500e3]

    hprimes = [72.0, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]

    altitude = 40e3:110e3
    density = [Vector{Float64}(undef, length(altitude)) for i = 1:length(segment_ranges)]
    collision_frequency = similar(density)
    for i in eachindex(segment_ranges)
        density[i] = waitprofile.(altitude, hprimes[i], betas[i])
        collision_frequency[i] = electroncollisionfrequency.(altitude)
    end

    b_mags = fill(50e-6, length(segment_ranges))
    b_dips = fill(π/2, length(segment_ranges))
    b_azs = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    frequency = 24e3

    output_ranges = collect(0:20e3:2000e3)

    input = TableInput()
    input.name = "table"
    input.description = "Test TableInput"
    input.datetime = Dates.now()

    input.segment_ranges = segment_ranges
    input.altitude = altitude
    input.density = density
    input.collision_frequency = collision_frequency
    input.b_mags= b_mags
    input.b_dips = b_dips
    input.b_azs = b_azs
    input.ground_sigmas = ground_sigmas
    input.ground_epsrs = ground_epsrs
    input.frequency = frequency
    input.output_ranges = output_ranges

    write && JSON.json("table.json", input; pretty=true)

    return input
end

function generate_table_power(power=1000.0; write=true)
    input = generate_table(write=false)
    input.power = power

    if power == 1000.0
        write && JSON.json("table.json", input; pretty=true)
    else
        write && JSON.json("table-power.json", input; pretty=true)
    end

    return input
end

function generate_batchtable()
    N = 2
    rep(v) = repeat(v, 1, N)

    # Waveguide definition
    nsegments = 4
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    b_mags = fill(50e-6, nsegments)
    b_dips = fill(π/2, nsegments)
    b_azs = fill(0.0, nsegments)
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    hprimes = [72.0, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]

    altitude = 40e3:110e3
    density = [Vector{Float64}(undef, length(altitude)) for i = 1:length(segment_ranges)]
    collision_frequency = similar(density)
    for i in eachindex(segment_ranges)
        density[i] = waitprofile.(altitude, hprimes[i], betas[i])
        collision_frequency[i] = electroncollisionfrequency.(altitude)
    end

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = collect(0:20e3:2000e3)

    binput = BatchInput{TableInput}()
    binput.name = "batchtable"
    binput.description = "Test BatchInput with TableInput"
    binput.datetime = Dates.now()

    inputs = Vector{TableInput}(undef, N)
    for i in eachindex(inputs)
        input = TableInput()
        input.name = "$i"
        input.description = "Test TableInput $i"
        input.datetime = Dates.now()

        input.segment_ranges = segment_ranges
        input.altitude = altitude
        input.density = density
        input.collision_frequency = collision_frequency
        input.b_mags= b_mags
        input.b_dips = b_dips
        input.b_azs = b_azs
        input.ground_sigmas = ground_sigmas
        input.ground_epsrs = ground_epsrs
        input.frequency = frequency
        input.output_ranges = output_ranges

        inputs[i] = input
    end

    binput.inputs = inputs

    JSON.json("batchtable.json", binput; pretty=true)

    return nothing
end

function generate_batchbasic(;write=true)
    N = 2
    rep(v) = repeat(v, 1, N)

    # Waveguide definition
    nsegments = 4
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    b_mags = fill(50e-6, nsegments)
    b_dips = fill(π/2, nsegments)
    b_azs = fill(0.0, nsegments)
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    hprimes = rep([72.0, 74, 75, 76])
    betas = rep([0.28, 0.30, 0.32, 0.35])

    frequency = 24e3

    output_ranges = collect(0:20e3:2000e3)

    binput = BatchInput{ExponentialInput}()
    binput.name = "batchbasic"
    binput.description = "Test BatchInput with ExponentialInput"
    binput.datetime = Dates.now()

    inputs = Vector{ExponentialInput}(undef, N)
    for i in eachindex(inputs)
        input = ExponentialInput()

        input.name = "$i"
        input.description = "ExponentialInput $i"
        input.datetime = binput.datetime
        input.segment_ranges = segment_ranges
        input.hprimes = hprimes[:,i]
        input.betas = betas[:,i]
        input.b_mags = b_mags
        input.b_dips = b_dips
        input.b_azs = b_azs
        input.ground_sigmas = ground_sigmas
        input.ground_epsrs = ground_epsrs
        input.frequency = frequency
        input.output_ranges = output_ranges

        inputs[i] = input
    end

    binput.inputs = inputs

    write && JSON.json("batchbasic.json", binput; pretty=true)

    return binput
end

function generate_batchbasic_power(power=1000.0; write=true)
    binput = generate_batchbasic(write=false)
    for input in binput.inputs
        input.power = power
    end

    if power == 1000.0
        write && JSON.json("batchbasic.json", binput; pretty=true)
    else
        write && JSON.json("batchbasic-power.json", binput; pretty=true)
    end

    return binput
end

read_basic() = JSON.parsefile("basic.json", ExponentialInput)
read_table() = JSON.parsefile("table.json", TableInput)
read_batchbasic() = JSON.parsefile("batchbasic.json", BatchInput{ExponentialInput})

# missing the ground_sigmas field
read_corrupted_basic() = JSON.parsefile("corrupted_basic.json", ExponentialInput)

function run_basic()
    s = LMP.parse("basic.json")
    output = LMP.buildrun(s)
    return output
end

function run_basic_fieldcomponent(fieldcomponent::Fields.Field)
    if fieldcomponent == Fields.Ez
        s = LMP.parse("basic.json")
    elseif fieldcomponent == Fields.E
        s = LMP.parse("basic-E.json")
    elseif fieldcomponent == Fields.EH
        s = LMP.parse("basic-EH.json")
    else
        error("Invalid field component: $fieldcomponent")
    end
    output = LMP.buildrun(s)
    return output
end

function run_basic_power(power::Float64)
    if power == 1000.0
        s = LMP.parse("basic.json")
    else
        s = LMP.parse("basic-power.json")
    end
    output = LMP.buildrun(s)
    return output
end

function run_table()
    s = LMP.parse("table.json")
    output = LMP.buildrun(s)
    return output
end

function run_table_power(power::Float64)
    if power == 1000.0
        s = LMP.parse("table.json")
    else
        s = LMP.parse("table-power.json")
    end
    output = LMP.buildrun(s)
    return output
end

function run_batchbasic()
    s = LMP.parse("batchbasic.json")
    output = LMP.buildrun(s)
    return output
end

function run_batchbasic_power(power::Float64)
    if power == 1000.0
        s = LMP.parse("batchbasic.json")
    else
        s = LMP.parse("batchbasic-power.json")
    end
    output = LMP.buildrun(s)
    return output
end

function run_batchbasicsave()
    s = LMP.parse("batchbasic.json")
    output = LMP.buildrunsave("batchbasictest.json", s)
    sres = LMP.parse("batchbasictest.json")
    return sres
end

function run_batchtablesave()
    s = LMP.parse("batchtable.json")
    output = LMP.buildrunsave("batchtabletest.json", s)
    sres = LMP.parse("batchtabletest.json")
    return sres
end

function basic_lmp()
    output = LMP.propagate("basic.json")

    # Compare to a `propagate` call with waveguide arguments
    tx = Transmitter(24e3)
    rx = GroundSampler(0:20e3:2000e3, Fields.Ez)

    # From `generate_basic()`
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    hprimes = [72.0, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]
    b_mags = fill(50e-6, length(segment_ranges))
    b_dips = fill(π/2, length(segment_ranges))
    b_azs = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    wvg = SegmentedWaveguide(
        [HomogeneousWaveguide(
            BField(b_mags[i], b_dips[i], b_azs[i]),
            Species(QE, ME, z->waitprofile(z, hprimes[i], betas[i]; cutoff_low=40e3), electroncollisionfrequency),
            Ground(ground_epsrs[i], ground_sigmas[i]),
            segment_ranges[i]
        ) for i = 1:4]
    )

    _, a, p = propagate(wvg, tx, rx)
    
    @test meanabsdiff(a, output.amplitude) < 0.1
    @test meanabsdiff(p, output.phase) < 0.01

    return output
end

@testset "IO.jl" begin
    @info "Testing IO"

    generate_basic()
    generate_basic_fieldcomponent(Fields.E)
    generate_basic_fieldcomponent(Fields.EH)
    generate_basic_power(500.0)
    generate_table()
    generate_table_power(500.0)
    generate_batchbasic()
    generate_batchbasic_power(500.0)
    generate_batchtable()

    @test LMP.iscomplete(read_basic())
    @test LMP.validlengths(read_basic())
    @test LMP.iscomplete(read_corrupted_basic()) == false
    @test LMP.parse("basic.json") isa ExponentialInput
    @test LMP.parse("basic.json", ExponentialInput) isa ExponentialInput
    @test LMP.parse("basic-E.json", ExponentialInput) isa ExponentialInput
    @test LMP.parse("basic-EH.json", ExponentialInput) isa ExponentialInput

    @test read_basic().fieldcomponent == Fields.Ez
    @test read_basic().power == 1000.0
    @test all(x->x.fieldcomponent == Fields.Ez, read_batchbasic().inputs)
    @test all(x->x.power == 1000.0, read_batchbasic().inputs)

    @test LMP.iscomplete(read_table())
    @test LMP.validlengths(read_table())
    @test LMP.parse("table.json") isa TableInput
    @test LMP.parse("table.json", TableInput) isa TableInput
    @test LMP.parse("table-power.json", TableInput) isa TableInput

    @test LMP.iscomplete(read_batchbasic())
    @test LMP.validlengths(read_batchbasic())
    @test LMP.parse("batchbasic.json") isa BatchInput{ExponentialInput}

    @info "  Running:"
    @info "    Segmented Wait ionospheres..."
    @test run_basic() isa BasicOutput
    @test run_basic_fieldcomponent(Fields.E) isa FieldsOutput
    @test run_basic_fieldcomponent(Fields.EH) isa FieldsOutput
    
    o1000 = run_basic_power(1000.0)
    o500 = run_basic_power(500.0)
    @test o1000 isa BasicOutput
    @test o500 isa BasicOutput
    @test all(>(2.9), o1000.amplitude[2:end] - o500.amplitude[2:end]) # ≈ 3 dB difference

    @info "    Segmented tabular ionospheres..."
    ot = run_table()
    ot500 = run_table_power(500.0)
    @test ot isa BasicOutput
    @test ot500 isa BasicOutput
    @test all(>(2.9), ot.amplitude[2:end] - ot500.amplitude[2:end]) # ≈ 3 dB difference
    @info "    Multiple segmented Wait ionospheres..."
    bo = run_batchbasic()
    bo500 = run_batchbasic_power(500.0)
    @test bo isa BatchOutput{BasicOutput}
    @test bo500 isa BatchOutput{BasicOutput}
    @test all(>(2.9), bo.outputs[1].amplitude[2:end] - bo500.outputs[1].amplitude[2:end]) # ≈ 3 dB difference
    @info "    Multiple segmented Wait ionospheres, appending..."
    @test run_batchbasicsave() isa BatchOutput{BasicOutput}
    @info "    Multiple segmented tabular ionospheres, appending."
    @test run_batchtablesave() isa BatchOutput{BasicOutput}

    @info "    `propagate` segmented Wait ionospheres"
    @test basic_lmp() isa BasicOutput

    isfile("basic.json") && rm("basic.json")
    isfile("basic-E.json") && rm("basic-E.json")
    isfile("basic-EH.json") && rm("basic-EH.json")
    isfile("basic-power.json") && rm("basic-power.json")
    isfile("table.json") && rm("table.json")
    isfile("table-power.json") && rm("table-power.json")
    isfile("batchbasic.json") && rm("batchbasic.json")
    isfile("batchbasic-power.json") && rm("batchbasic-power.json")
    isfile("batchtable.json") && rm("batchtable.json")
    isfile("batchbasictest.json") && rm("batchbasictest.json")
    isfile("batchtabletest.json") && rm("batchtabletest.json")
    isfile("basic_output.json") && rm("basic_output.json")
end
