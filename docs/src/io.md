# File-based I/O

LongwaveModePropagator provides basic propagation capabilities
interfacing through [JSON](https://www.json.org/json-en.html) files.
Julia is still needed to run the model, but scenarios can otherwise
be defined and analyzed from e.g. Matlab or Python.

Examples of writing and reading compatible JSON files are provided below
for [Matlab](@ref matlab_json) and [Python](@ref python_json).

Let's load the necessary packages.

```@example io
using Dates
using JSON3

using LongwaveModePropagator
```

Throughout the examples, we'll also define `LMP` as shorthand for
`LongwaveModePropagator` when accessing functions from the package that
aren't exported.

```@example io
const LMP = LongwaveModePropagator
```

## Inputs

There are two primary ways to define [`LongwaveModePropagator.Input`](@ref)s,
the abstract supertype for inputing information to the model.

### [ExponentialInput](@id basicinput_io)
 
The first is known as a [`ExponentialInput`](@ref). It defines the ionosphere using
Wait and Spies ``h'`` and ``\beta`` parameters.

```@docs; canonical=false
ExponentialInput
```

The fields that are vectors allow the definition of a [`SegmentedWaveguide`](@ref)
where each element of the vector is its own [`HomogeneousWaveguide`](@ref) segment.
Single element vectors are treated as a single `HomogeneousWaveguide`.

To show the equivalent JSON format, we'll build a simple, homogeneous ionosphere
`ExponentialInput`. It's defined as a `mutable struct`, so it is simple to specify the
fields one by one.

```@example io
input = ExponentialInput()
input.name = "basic"
input.description = "Test ExponentialInput"
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
```

Here it is formatted as JSON.

```@example io
json_str = JSON3.pretty(input)
```

Let's also save this to a file.

```@example io
json_str = JSON3.write(input)

filename = "basic.json"

open(filename,"w") do f
    write(f, json_str)
end
```

### [TableInput](@id tableinput_io)
 
The second input is the [`TableInput`](@ref). This type defines the ionosphere using
a tabular input of number density and collision frequency as a function of altitude.
These tables are then cubic spline interpolated when integrating the ionosphere
reflection coefficient and wavefields. See also [interpolating functions](@ref interpolating-functions).

```@docs; canonical=false
TableInput
```

Again we'll construct a `TableInput` to look at the JSON

```@example io
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
```

## [BatchInput](@id batchinput_io)
 
Both the `ExponentialInput` and `TableInput` types can be collected together in a
[`BatchInput`](@ref) which has fields
for a `name`, `description`, `datetime`, and vector of `Inputs`.
This is useful for keeping a set of scenarios together.
See [test/IO.jl](https://github.com/fgasdia/LongwaveModePropagator.jl/test/IO.jl) for additional help on how these should be formatted.
 
## Running the model from a JSON file

To run the model, [`propagate`](@ref) accepts a filename input
(see the help for optional arguments). However, instead of returning a
tuple of complex electric field, amplitude, and phase, it returns an `Output`
type. Additionally, it saves the output to a JSON file.

```@example io
output = propagate(filename);
```

## Outputs

There are only two [`LongwaveModePropagator.Output`](@ref) types:
[`BasicOutput`](@ref) and [`BatchOutput`](@ref).
Both `ExponentialInput`s and `TableInput`s create `BasicOutput`s, but the `BatchInput`
creates a `BatchOutput`.

```@docs; canonical=false
BasicOutput
```

The `name` and `description` fields were copied from the input. `datetime`
is when the model was run.

```@example io
output
```

Not surprisingly, a `BatchOutput` is simply a container holding a `Vector` of
`BasicOutput`s, as well as some additional metadata from the corresponding `BatchInput`.

## [JSON I/O from Matlab](@id matlab_json)

Here's an example of how to encode the above `ExponentialInput` to JSON and
decode the output using [Matlab](https://www.mathworks.com/help/matlab/json-format.html).
It's also in the file [io.m](https://github.com/fgasdia/LongwaveModePropagator.jl/docs/src/io.m).

```matlab
% Matlab script

input.name = "basic";
input.description = "Test ExponentialInput";
input.datetime = '2020-12-28T21:06:50.501';

input.segment_ranges = {0.0};
input.hprimes = {75};
input.betas = {0.35};
input.b_mags = {50e-6};
input.b_dips = {pi/2};
input.b_azs = {0.0};
input.ground_sigmas = {0.001};
input.ground_epsrs = {4};
input.frequency = 24e3;
input.output_ranges = 0:100e3:1000e3;

json_str = jsonencode(input);

fid = fopen('basic_matlab.json', 'w');
fwrite(fid, json_str, 'char');
fclose(fid);
```

Matlab was used to generate the file
[`basic_matlab.json`](https://github.com/fgasdia/LongwaveModePropagator.jl/docs/src/basic_matlab.json).
We can confirm it's parsed correctly by using the internal
LongwaveModePropagator function [`LongwaveModePropagator.parse`](@ref), which attempts to parse
JSON files into recognized input and output formats.

```@example io
matlab_input = LMP.parse("basic_matlab.json")
```

Let's run it.

```@example io
matlab_output = propagate("basic_matlab.json");
```

!!! note
    It is possible to run Julia as a script, calling it directly from a terminal
    (see [the docs](https://docs.julialang.org/en/v1/manual/getting-started/)), but
    it is probably easiest to just run the code from the REPL.

    1. `cd` in a terminal to your working directory
    2. Start up Julia: `julia`
    3. It's recommended to `] activate .`, generating a new environment in this directory
    4. If necessary, install LongwaveModePropagator.jl: `] add LongwaveModePropagator`
    5. `using LongwaveModePropagator`
    6. `propagate("basic_matlab.json")`

    If this has been done before in the same directory, then just do steps 1, 2, 5, 6.

Reading the results file from Matlab is relatively simple.

```matlab
% Matlab script

filename = 'basic_matlab_output.json';

text = fileread(filename);

output = jsondecode(text);
```

## [JSON I/O from Python](@id python_json)

Here is similar code for [Python](https://docs.python.org/3/library/json.html),
also available in the file [io.py](https://github.com/fgasdia/LongwaveModePropagator.jl/docs/src/io.py).

```python
# Python script

import json
import datetime
import numpy as np

input = dict()

input['name'] = "basic"
input['description'] = "Test ExponentialInput"
input['datetime'] = datetime.datetime.now().isoformat()[:-3]

input['segment_ranges'] = [0.0]
input['hprimes'] = [75]
input['betas'] = [0.35]
input['b_mags'] = [50e-6]
input['b_dips'] = [np.pi/2]
input['b_azs'] = [0.0]
input['ground_sigmas'] = [0.001]
input['ground_epsrs'] = [4]
input['frequency'] = 24e3
input['output_ranges'] = np.arange(0, 1000e3, 100e3).tolist()

json_str = json.dumps(input)

with open('basic_python.json', 'w') as file:
    file.write(json_str)
```

Let's ensure that the JSON file is correctly parsed.

```@example io
python_input = LMP.parse("basic_python.json")
```

And run the file.

```@example io
python_output = propagate("basic_python.json");
```

To read the results:

```python
with open('basic_python_output.json', 'r') as file:
    output = json.load(file)
```
