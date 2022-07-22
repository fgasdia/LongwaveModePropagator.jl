# Public Interface

Documentation of `LongwaveModePropagator.jl`'s exported structs and functions.

**Contents**

- [Basic Functions](@ref)
- [Mode Finder](@ref)
- [Geophysics](@ref)
- [Samplers](@ref)
- [Emitters](@ref)
- [Waveguides](@ref)
- [IO](@ref)

### Basic Functions

```@docs
propagate
LMPParams
```

### Mode Finder

```@docs
findmodes
PhysicalModeEquation
setea
IntegrationParams
```

### EigenAngle

```@docs
attenuation
phasevelocity
referencetoground
```

### Geophysics

```@docs
BField
Species
Ground
GROUND
waitprofile
electroncollisionfrequency
ioncollisionfrequency
```

### Samplers

```@docs
Sampler
GroundSampler
Fields
```

### Emitters

```@docs
Transmitter
Dipole
VerticalDipole
HorizontalDipole
inclination
azimuth
```

### Waveguides

```@docs
HomogeneousWaveguide
SegmentedWaveguide
```

### IO

```@docs
BasicInput
TableInput
BatchInput
BasicOutput
BatchOutput
```
