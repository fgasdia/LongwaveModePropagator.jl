% Matlab script to create a JSON file formated to be parsed
% by LongwaveModePropagator.jl as a `BasicInput`.
% 
% See examples/io.jl.

input.name = "basic";
input.description = "Test BasicInput";
input.datetime = '2020-12-28T21:06:50.501';

input.segment_ranges = [0.0];
input.hprimes = [75];
input.betas = [0.35];
input.b_mags = [50e-6];
input.b_dips = [pi/2];
input.b_azs = [0.0];
input.ground_sigmas = [0.001];
input.ground_epsrs = [4];
input.frequency = 24e3;
input.output_ranges = 0:100e3:1000e3;

json_str = jsonencode(input);

fid = fopen('basic_matlab.json', 'w');
fwrite(fid, json_str, 'char');
fclose(fid);