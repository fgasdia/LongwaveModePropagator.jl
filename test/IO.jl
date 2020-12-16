function generate_basic()
    # Waveguide definition
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    hprimes = [72.0, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]
    b_mags = fill(50e-6, length(segment_ranges))
    b_dips = fill(π/2, length(segment_ranges))
    b_azs = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = collect(0:20e3:2000e3)

    input = BasicInput()
    input.name = "basic"
    input.description = "Test BasicInput"
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

    json_str = JSON3.write(input)

    open("basic.json","w") do f
        write(f, json_str)
    end

    return nothing
end

function generate_table()
    # Waveguide definition
    segment_ranges = [0, 500e3, 1000e3, 1500e3]

    hprimes = [72.0, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]

    altitude = 40e3:110e3
    density = [Vector{Float64}(undef, length(altitude)) for i = 1:length(segment_ranges)]
    collision_frequency = similar(density)
    for i in eachindex(segment_ranges)
        density[i] = LMP.waitprofile.(altitude, hprimes[i], betas[i])
        collision_frequency[i] = LMP.electroncollisionfrequency.(altitude)
    end

    b_mags = fill(50e-6, length(segment_ranges))
    b_dips = fill(π/2, length(segment_ranges))
    b_azs = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    # Transmitter
    frequency = 24e3

    # Outputs
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

    json_str = JSON3.write(input)

    open("table.json","w") do f
        write(f, json_str)
    end

    return nothing
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
        density[i] = LMP.waitprofile.(altitude, hprimes[i], betas[i])
        collision_frequency[i] = LMP.electroncollisionfrequency.(altitude)
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

    json_str = JSON3.write(binput)

    open("batchtable.json","w") do f
        write(f, json_str)
    end

    return nothing
end

function generate_batchbasic()
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

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = collect(0:20e3:2000e3)

    binput = BatchInput{BasicInput}()
    binput.name = "batchbasic"
    binput.description = "Test BatchInput with BasicInput"
    binput.datetime = Dates.now()

    inputs = Vector{BasicInput}(undef, N)
    for i in eachindex(inputs)
        input = BasicInput()

        input.name = "$i"
        input.description = "BasicInput $i"
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

    json_str = JSON3.write(binput)

    open("batchbasic.json","w") do f
        write(f, json_str)
    end

    return nothing
end

function read_basic()
    open("basic.json","r") do f
        v = JSON3.read(f, BasicInput)
        return v
    end
end

function read_table()
    open("table.json", "r") do f
        v = JSON3.read(f, TableInput)
        return v
    end
end

function read_batchbasic()
    open("batchbasic.json","r") do f
        v = JSON3.read(f, BatchInput{BasicInput})
        return v
    end
end

function read_corrupted_basic()
    # missing the ground_sigmas field
    open("corrupted_basic.json","r") do f
        v = JSON3.read(f, BasicInput)
        return v
    end
end

function run_basic()
    s = LMP.parse("basic.json")
    output = LMP.buildrun(s)
    return output
end

function run_table()
    s = LMP.parse("table.json")
    output = LMP.buildrun(s)
    return output
end

function run_batchbasic()
    s = LMP.parse("batchbasic.json")
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
    return output
end

@testset "IO.jl" begin
    @info "Testing IO"

    generate_basic()
    generate_table()
    generate_batchbasic()
    generate_batchtable()

    @test LMP.iscomplete(read_basic())
    @test LMP.validlengths(read_basic())
    @test LMP.iscomplete(read_corrupted_basic()) == false
    @test LMP.parse("basic.json") isa BasicInput
    @test LMP.parse("basic.json", BasicInput) isa BasicInput

    @test LMP.iscomplete(read_table())
    @test LMP.validlengths(read_table())
    @test LMP.parse("table.json") isa TableInput
    @test LMP.parse("table.json", TableInput) isa TableInput

    @test LMP.iscomplete(read_batchbasic())
    @test LMP.validlengths(read_batchbasic())
    @test LMP.parse("batchbasic.json") isa BatchInput{BasicInput}

    @info "  Running:"
    @info "    Segmented Wait ionospheres..."
    @test run_basic() isa BasicOutput
    @info "    Segmented tabular ionospheres..."
    @test run_table() isa BasicOutput
    @info "    Multiple segmented Wait ionospheres..."
    @test run_batchbasic() isa BatchOutput{BasicOutput}
    @info "    Multiple segmented Wait ionospheres, appending..."
    @test run_batchbasicsave() isa BatchOutput{BasicOutput}
    @info "    Multiple segmented tabular ionospheres, appending."
    @test run_batchtablesave() isa BatchOutput{BasicOutput}

    @info "    `propagate` segmented Wait ionospheres"
    @test basic_lmp() isa BasicOutput

    isfile("basic.json") && rm("basic.json")
    isfile("table.json") && rm("table.json")
    isfile("batchbasic.json") && rm("batchbasic.json")
    isfile("batchtable.json") && rm("batchtable.json")
    isfile("batchbasictest.json") && rm("batchbasictest.json")
    isfile("batchtabletest.json") && rm("batchtabletest.json")
    isfile("basic_output.json") && rm("basic_output.json")
end
