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
    s = LWMS.parse("basic.json")
    output = LWMS.buildrun(s)
    return output
end

function run_batchbasic()
    s = LWMS.parse("batchbasic.json")
    output = LWMS.buildrun(s)
    return output
end

function run_batchbasicsave()
    s = LWMS.parse("batchbasic.json")
    output = LWMS.buildrunsave("batchbasictest.json", s)
    sres = LWMS.parse("batchbasictest.json")
    return sres
end

function test_bpm()
    E, amp, phase = LWMS.bpm("basic.json")
end

@testset "IO.jl" begin
    @info "Testing IO"

    generate_basic()
    generate_batchbasic()

    @test LWMS.iscomplete(read_basic())
    @test LWMS.validlengths(read_basic())
    @test LWMS.iscomplete(read_corrupted_basic()) == false
    @test LWMS.parse("basic.json") isa BasicInput
    @test LWMS.parse("basic.json", BasicInput) isa BasicInput

    @test LWMS.iscomplete(read_batchbasic())
    @test LWMS.validlengths(read_batchbasic())
    @test LWMS.parse("batchbasic.json") isa BatchInput{BasicInput}

    @test run_basic() isa BasicOutput
    @test run_batchbasic() isa BatchOutput{BasicOutput}
    @test run_batchbasicsave() isa BatchOutput{BasicOutput}

    @test_skip test_bpm()

    isfile("basic.json") && rm("basic.json")
    isfile("batchbasic.json") && rm("batchbasic.json")
end
