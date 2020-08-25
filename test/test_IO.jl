function generate_basic()
    # Waveguide definition
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    hprimes = [72, 74, 75, 76]
    betas = [0.28, 0.30, 0.32, 0.35]
    b_mag = fill(50e-6, length(segment_ranges))
    b_dip = fill(90, length(segment_ranges))
    b_az = fill(0, length(segment_ranges))
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsr = [4, 4, 10, 10]

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = 0:20e3:2000e3

    json_str = JSON3.write(BasicInput(segment_ranges, hprimes, betas,
        b_mag, b_dip, b_az,
        ground_sigmas, ground_epsr, frequency, output_ranges))

    open(joinpath("test","basic.json"),"w") do f
        write(f, json_str)
    end
end

function read_basic()
    open(joinpath("test", "basic.json"),"r") do f
        v = JSON3.read(f, BasicInput)
        return v
    end
end

function read_corrupted_basic()
    open(joinpath("test", "corrupted_basic.json"),"r") do f
        v = JSON3.read(f, BasicInput)
        return v
    end
end

@testset "Testing IO" begin
    @test iscomplete(read_basic()) == true
    @test isvalid(read_basic()) == true
    @test iscomplete(read_corrupted_basic()) == false
end
