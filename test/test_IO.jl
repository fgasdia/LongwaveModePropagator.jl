function generate_basic()
    # Waveguide definition
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    hprimes = [72.0, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]
    b_mag = fill(50e-6, length(segment_ranges))  # TODO: rename to be consistent (plural?)
    b_dip = fill(90.0, length(segment_ranges))
    b_az = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsr = [4, 4, 10, 10]

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
    input.b_mag = b_mag
    input.b_dip = b_dip
    input.b_az = b_az
    input.ground_sigmas = ground_sigmas
    input.ground_epsr = ground_epsr
    input.frequency = frequency
    input.output_ranges = output_ranges

    json_str = JSON3.write(input)

    open("basic.json","w") do f
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

function read_corrupted_basic()
    # missing the ground_sigmas field
    open("corrupted_basic.json","r") do f
        v = JSON3.read(f, BasicInput)
        return v
    end
end

function test_bpm()
    E, amp, phase = LWMS.bpm("basic.json")
end

@testset "Testing IO" begin
    generate_basic()

    @test LWMS.iscomplete(read_basic()) == true
    @test LWMS.validlengths(read_basic()) == true
    @test LWMS.iscomplete(read_corrupted_basic()) == false

    @test_skip test_bpm()

    isfile("basic.json") && rm("basic.json")
end
