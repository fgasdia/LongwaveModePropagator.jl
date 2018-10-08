"""
Main program/scratch file
"""
module LWMS

using Dates
using Proj4

export earthradius
export Inputs, LWMSModes

const earthradius = 6369.427  # km, because it is used with heights

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

struct LWMSModes
    ranger::Array{Float64}  # (undef, 2)
    rangei::Array{Float64}  # (undef, 2)
    atnmax::Float64
    deltathetathreshold::Float64
    height::Float64
    eigen::Array{ComplexF64}
    numeigen::Integer
end

function main()
    program_start_time = Dates.now()
    print("Program started: $program_start_time")

    # Define test parameters
    xmtr = Transmitter("NAA", 44.646394, -67.281069, 0.0, 24.0)
    rcvr = Receiver("Boulder", 40.014984, -105.270546, 0.0)

    # Get bearing angles
    # TODO: Make bearing on (0, 360)?
    # TODO: Confirm with LWPC that bearing is from TX to RX
    earthprojection = Projection("+proj=longlat +a=6366200 +b=6366200")
    wgs84 = Projection("+proj=longlat +datum=WGS84 +no_defs")
    distance, bearing = Proj4._geod_inverse(Proj4._geod(earthprojection),
        [xmtr.lon, xmtr.lat], [rcvr.lon, rcvr.lat])[1:2]

    # Round to nearest 100 km
    maxrange = round(distance, digits=-5, RoundUp)

    # Get modes

end

end
