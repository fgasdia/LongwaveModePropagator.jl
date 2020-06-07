"""
this function takes just 1 mode conversion step

"""
function modeconversion(previous_wavefields::Wavefields{T},
                        wavefields::Wavefields{T}, adjwavefields::Wavefields{T}) where {T}

    @assert numheights(previous_wavefields) == numheights(wavefields)
    zs = heights(wavefields)

    # TODO: Assuming `length(zs)` is always the same, we can reuse `product`
    product = Vector{ComplexF64}(undef, numheights(wavefields))

    modes = eigenangles(wavefields)
    adjmodes = eigenangles(adjwavefields)
    prevmodes = eigenangles(previous_wavefields)

    nummodes = length(modes)
    numadjmodes = length(adjmodes)
    numprevmodes = length(prevmodes)
    @assert nummodes == numadjmodes  # for correct physics, need modes == adjmodes


    # Calculate normalization terms
    N = Vector{ComplexF64}(undef, nummodes)
    for m in eachindex(modes)
        for i in eachindex(zs)
            @inbounds f = wavefields[m][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
            @inbounds g = SVector{4,ComplexF64}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
            product[i] = adjoint(g)*f
        end

        N[m] = integrate(zs, product, RombergEven())
        # @test isapprox(N[m], trapz(zs, product), rtol=1e-3)
        # @test isapprox(N[m], romberg(zs, product), rtol=1e-3)
    end

    I = Matrix{ComplexF64}(undef, numadjmodes, numprevmodes)
    for k in eachindex(prevmodes)
        for m in eachindex(adjmodes)
            for i in eachindex(zs)
                @inbounds f = previous_wavefields[k][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
                @inbounds g = SVector{4,ComplexF64}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
                product[i] = adjoint(g)*f
            end
            I[m,k] = integrate(zs, product, RombergEven())
        end
    end

    # Total conversion
    # a = Matrix{ComplexF64}(undef, numprevmodes, numadjmodes)
    # for k in eachindex(prevmodes)
    #     for m in eachindex(adjmodes)
    #         a[k,m] = I[m,k]/N[m]
    #     end
    # end
    a = I./N

    return a
end
