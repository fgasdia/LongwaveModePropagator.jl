using Dates
using Proj4

using LongwaveModeSolver

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
    modesolver()
end
