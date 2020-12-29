# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 11:14:54 2020

Python script to create a JSON file formated to be parsed by
LongwaveModePropagator.jl as a `BasicInput`.

See examples/io.jl.

@author: forrest
"""

import json
import datetime
import numpy as np

# We'll use a dictionary for simplicity, but look at JSON docs to build custom
# JSONEncoders with a class for the Input types in LongwaveModePropagator.jl.

input = dict()

input['name'] = "basic"
input['description'] = "Test BasicInput"
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
