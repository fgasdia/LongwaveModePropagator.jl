"""
    readlog(file)

Read standard LWPC `.log` files and return `distance`, `amplitude`, and `phase` vectors.

The amplitude and phase at the receiver is not returned by this function.
"""
function readlog(file)
    # Search for first and last line of data
    lines = readlines(file)
    firstdataline = findfirst(startswith.(lines, "  dist   amplitude  phase")) + 1
    lastdataline = findfirst(startswith.(lines, "nc nrpt bearng")) - 1
    numdatalines = lastdataline - firstdataline + 1

    # Read log file
    raw = CSV.File(file; datarow=firstdataline, limit=numdatalines,
                   delim=' ', ignorerepeated=true, header=false,
                   silencewarnings=true, threaded=false)

    # Stack the columns together
    dist = vcat(raw.Column1, raw.Column4, raw.Column7)
    amplitude = vcat(raw.Column2, raw.Column5, raw.Column8)
    phase = vcat(raw.Column3, raw.Column6, raw.Column9)

    # Clean up end of the last column
    delidxs = ismissing.(dist)  # assuming these are at the end
    delidxs[findfirst(delidxs)-1] = true  # last valid entry is at receiver distance
    deleteat!(dist, delidxs)
    deleteat!(amplitude, delidxs)
    deleteat!(phase, delidxs)

    # If phase gets above 9999 deg in log file, there is no space between amplitude and phase
    if count(ismissing.(amplitude)) != count(ismissing.(phase))
        for i in eachindex(phase)
            if ismissing(phase[i])
                phase[i] = parse(Float64, amplitude[i][end-9:end])  # phase is 10000.0000
                amplitude[i] = parse(Float64, amplitude[i][1:end-10])
            end
        end
        # Other elements in the same column will also be string type
        for i in eachindex(amplitude)
            if amplitude[i] isa String
                amplitude[i] = parse(Float64, amplitude[i])
            end
        end
    end

    return dist, amplitude, phase
end
