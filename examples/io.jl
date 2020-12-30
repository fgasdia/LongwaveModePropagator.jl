# # File-based I/O
#
# LongwaveModePropagator.jl provides basic propagation capabilities
# through input/output JSON files. Julia is still needed to run the
# model, but scenarios can otherwise be defined and analyzed from
# e.g. Matlab or Python.
#
# Examples of writing and reading compatible JSON files are provided below
# for [Matlab](@ref matlab_json) and [Python](@ref python_json).
#
# Let's load the necessary packages.

using Dates
using JSON3
using Plots
using DisplayAs  #hide

using ..LongwaveModePropagator
const LMP = LongwaveModePropagator
nothing  #hide

# ## Inputs
#
# There are two primary ways to define `Input`s.
#
# The first is known as a [`BasicInput`](@ref). It defines the ionosphere using
# Wait and Spies ``h'`` and ``\beta`` parameters.
#
# It contains the fields
#
# - `name::String`
# - `description::String`
# - `datetime::DateTime`
# - `segment_ranges::Vector{Float64}`: distance from transmitter to the beginning of each `HomogeneousWaveguide` segment in meters.
# - `hprimes::Vector{Float64}`: Wait's ``h'`` parameter for each `HomogeneousWaveguide` segment.
# - `betas::Vector{Float64}`: Wait's ``\beta`` parameter for each `HomogeneousWaveguide` segment.
# - `b_mags::Vector{Float64}`: magnetic field magnitude for each `HomogeneousWaveguide` segment.
# - `b_dips::Vector{Float64}`: magnetic field dip angles in radians for each `HomogeneousWaveguide` segment.
# - `b_azs::Vector{Float64}`: magnetic field azimuth in radians "east" of the propagation direction for each `HomogeneousWaveguide` segment.
# - `ground_sigmas::Vector{Float64}`: ground conductivity in Siemens per meter for each `HomogeneousWaveguide` segment.
# - `ground_epsrs::Vector{Int}`: ground relative permittivity for each `HomogeneousWaveguide` segment.
# - `frequency::Float64`: transmitter frequency in Hertz.
# - `output_ranges::Vector{Float64}`: distances from the transmitter at which the field will be calculated.
#
# The fields that are vectors allow the definition of a `SegmentedWaveguide`.
#
# To show the equivalent JSON format, we'll build a simple, homoegenous ionosphere
# `BasicInput`. It's defined as a `mutable struct`, so it is simple to specify the
# fields one by one.

input = BasicInput()
input.name = "basic"
input.description = "Test BasicInput"
input.datetime = DateTime("2020-12-29T05:00:00.000")  # usually `Dates.now()`

input.segment_ranges = [0.0]
input.hprimes = [75]
input.betas = [0.35]
input.b_mags= [50e-6]
input.b_dips = [π/2]
input.b_azs = [0.0]
input.ground_sigmas = [0.001]
input.ground_epsrs = [4]
input.frequency = 24e3
input.output_ranges = collect(0:100e3:1000e3)

json_str = JSON3.pretty(input)

# Let's also save this to a file.

json_str = JSON3.write(input)

examples_dir = joinpath("..", "..", "..", "examples")
filename = joinpath(examples_dir, "basic.json")

open(filename,"w") do f
    write(f, json_str)
end
nothing  #hide

# The second input is the `TableInput`. This type defines the ionosphere using
# a tabular input of number density and collision frequency over altitude.
# These tables are then _linearly_ interpolated when integrating the ionosphere
# reflection coefficient and wavefields (so it's better if these tables are fairly
# dense).
#
# The fields of the `TableInput` are
#
# - `name::String`
# - `description::String`
# - `datetime::DateTime`
# - `segment_ranges::Vector{Float64}`: distance from transmitter to the beginning of each `HomogeneousWaveguide` segment in meters.
# - `altitude::Vector{Float64}`: altitude above ground in meters for which the `density` and `collision_frequency` profiles are specified.
# - `density::Vector{Float64}`: electron density at each `altitude` in ``m⁻³``.
# - `collision_frequency::Vector{Float64}`: electron-ion collision frequency at each `altitude` in ``s⁻¹``.
# - `b_dips::Vector{Float64}`: magnetic field dip angles in radians for each `HomogeneousWaveguide` segment.
# - `b_azs::Vector{Float64}`: magnetic field azimuth in radians "east" of the propagation direction for each `HomogeneousWaveguide` segment.
# - `ground_sigmas::Vector{Float64}`: ground conductivity in Siemens per meter for each `HomogeneousWaveguide` segment.
# - `ground_epsrs::Vector{Int}`: ground relative permittivity for each `HomogeneousWaveguide` segment.
# - `frequency::Float64`: transmitter frequency in Hertz.
# - `output_ranges::Vector{Float64}`: distances from the transmitter at which the field will be calculated.
#
# Again we'll construct a `TableInput` to look at the JSON

tinput = TableInput()
tinput.name = "table"
tinput.description = "Test TableInput"
tinput.datetime = DateTime("2020-12-29T05:00:00.000")

tinput.segment_ranges = [0.0]
tinput.altitude = collect(50e3:5e3:100e3)
tinput.density = [waitprofile.(tinput.altitude, 75, 0.3)]
tinput.collision_frequency = [electroncollisionfrequency.(tinput.altitude)]
tinput.b_mags= [50e-6]
tinput.b_dips = [π/2]
tinput.b_azs = [0.0]
tinput.ground_sigmas = [0.001]
tinput.ground_epsrs = [4]
tinput.frequency = 24e3
tinput.output_ranges = collect(0:100e3:1000e3)

json_str = JSON3.pretty(tinput)

# Both of these input types can be collected together in a `BatchInput` which has fields
# for a `name`, `description`, `datetime`, and vector of inputs.
# This is useful for keeping a set of scenarios together.
# See the "test/IO.jl" file for additional help on how these should be formatted.
#
# To run the model, [`propagate`](@ref) accepts a filename input
# (see the help for optional arguments), however, instead of returning a
# tuple of complex electric field, amplitude, and phase, it returns an output
# type. It also saves the output to a JSON file.

output = propagate(filename);

# ## Outputs
#
# There are only two `Output` types: [`BasicOutput`](@ref) and [`BatchOutput`](@ref).
# Both `BasicInput` and `TableInput`s create `BasicOutput`s, but the `BatchInput`
# creates a `BatchOutput`.
#
# The `BasicOutput` contains fields for
#
# - `name::String`
# - `description::String`
# - `datetime::DateTime`
# - `output_ranges::Vector{Float64}`
# - `amplitude::Vector{Float64}`
# - `phase::Vector{Float64}`
#
# where `name` and `description` are directly copied from the input and `datetime`
# is when the model was run.

output

# ## [JSON I/O from Matlab](@id matlab_json)
#
# Here's an example of how to encode the above `BasicInput` to JSON and
# decode the output using Matlab.
# It's also in the file [io.m](@__REPO_ROOT_URL__/examples/io.m).
#
# ```matlab
# % Matlab script
#
# input.name = "basic";
# input.description = "Test BasicInput";
# input.datetime = '2020-12-28T21:06:50.501';
#
# input.segment_ranges = {0.0};
# input.hprimes = {75};
# input.betas = {0.35};
# input.b_mags = {50e-6};
# input.b_dips = {pi/2};
# input.b_azs = {0.0};
# input.ground_sigmas = {0.001};
# input.ground_epsrs = {4};
# input.frequency = 24e3;
# input.output_ranges = 0:100e3:1000e3;
#
# json_str = jsonencode(input);
#
# fid = fopen('basic_matlab.json', 'w');
# fwrite(fid, json_str, 'char');
# fclose(fid);
# ```
#
# Matlab was used to generate the file `basic_matlab.json`.
# We can confirm it's parsed correctly by using the internal
# LongwaveModePropagator function [`LongwaveModePropagator.parse`](@ref), which attempts to parse
# JSON files into recognized input and output formats.

matlab_input = LMP.parse(joinpath(examples_dir, "basic_matlab.json"))

# Let's run it.

matlab_output = propagate(joinpath(examples_dir, "basic_matlab.json"));

# !!! note
#     It is possible to run Julia as a script, calling it directly from a terminal
#     (see [the docs](https://docs.julialang.org/en/v1/manual/getting-started/)), but
#     it is probably easiest to just run the code from the REPL.
#
#     1. `cd` in a terminal to your working directory
#     2. Start up Julia: `julia`
#     3. It's recommended to `] activate .`, generating a new environment in this directory
#     4. If necessary, install LongwaveModePropagator.jl: `] add LongwaveModePropagator`
#     5. `using LongwaveModePropagator`
#     6. `propagate("basic_matlab.json")`
#
#     If this has been done before in the same directory, then just do steps 1, 2, 5, 6.
#
# Reading the results file from Matlab is relatively simple.
#
# ```matlab
# filename = 'basic_matlab_output.json';
#
# text = fileread(filename);
#
# output = jsondecode(text);
# ```
#
# ## [JSON I/O from Python](@id python_json)
#
# Here is similar code for Python, also available in the file
# [io.py](@__REPO_ROOT_URL__/examples/io.py).
#
# ```python
# # Python script
#
# import json
# import datetime
# import numpy as np
#
# input = dict()
#
# input['name'] = "basic"
# input['description'] = "Test BasicInput"
# input['datetime'] = datetime.datetime.now().isoformat()[:-3]
#
# input['segment_ranges'] = [0.0]
# input['hprimes'] = [75]
# input['betas'] = [0.35]
# input['b_mags'] = [50e-6]
# input['b_dips'] = [np.pi/2]
# input['b_azs'] = [0.0]
# input['ground_sigmas'] = [0.001]
# input['ground_epsrs'] = [4]
# input['frequency'] = 24e3
# input['output_ranges'] = np.arange(0, 1000e3, 100e3).tolist()
#
# json_str = json.dumps(input)
#
# with open('basic_python.json', 'w') as file:
#     file.write(json_str)
# ```
#
# Let's ensure that the JSON file is correctly parsed.

python_input = LMP.parse(joinpath(examples_dir, "basic_python.json"))

# And run the file.

python_output = propagate(joinpath(examples_dir, "basic_python.json"));

# To read the results:
#
# ```python
# with open('basic_python_output.json', 'r') as file:
#     output = json.load(file)
# ```
