export Receiver, Transmitter, Inputs

struct Receiver
    name::String
    lat::Float64
    lon::Float64
    alt::Float64
end

struct Transmitter
    name::String
    lat::Float64
    lon::Float64
    alt::Float64
    freq::Float64
end

# TEMP: Only need to be mutable and have internal constructor for incomplete initialization
mutable struct Inputs
    transmitter::Transmitter
    receiver::Receiver
    maxrange::Float64
    deltarange::Float64
    topheight::Float64
    bottomheight::Float64
    Inputs() = new()
end
