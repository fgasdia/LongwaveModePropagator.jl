"""
this function takes just 1 mode conversion step

returns a matrix...
"""
function modeconversion(previous_wavefields::Wavefields{T,T2},
                        wavefields::Wavefields{T,T2}, adjwavefields::Wavefields{T,T2},
                        transmitter_slab=false) where {T,T2}

    @assert numheights(previous_wavefields) == numheights(wavefields)
    zs = heights(wavefields)

    # TODO: Assuming `length(zs)` is always the same, we can reuse `product`
    product = Vector{T}(undef, numheights(wavefields))

    # Calculate normalization terms
    N = Vector{T}(undef, numeigenangles(wavefields))
    for m in eachindex(eigenangles(wavefields))
        for i in eachindex(zs)
            @inbounds f = wavefields[m][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
            @inbounds g = SVector{4,T}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
            product[i] = adjoint(g)*f
        end

        N[m] = integrate(zs, product, RombergEven())
        # @test isapprox(N[m], trapz(zs, product), rtol=1e-3)
        # @test isapprox(N[m], romberg(zs, product), rtol=1e-3)
    end

    # if transmitter_slab
    #     a = zeros(T, length(eas), length(eas))
    #     for i in diagind(a)
    #         @inbounds a[i] = one(ComplexF64)
    #     end
    #
    #     # TODO: a = I ?  # UniformScaling identity matrix, probably need to pull this out for type consistency
    # else
        I = Matrix{T}(undef, numeigenangles(adjwavefields), numeigenangles(previous_wavefields))
        for k in eachindex(eigenangles(previous_wavefields))
            for m in eachindex(eigenangles(adjwavefields))
                for i in eachindex(heights(adjwavefields))
                    @inbounds f = previous_wavefields[k][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
                    @inbounds g = SVector{4,T}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
                    product[i] = adjoint(g)*f
                end
                I[m,k] = integrate(zs, product, RombergEven())
            end
        end

        # Total conversion
        for k in eachindex(eigenangles(previous_wavefields))
            for m in eachindex(eigenangles(adjwavefields))
                a[k,m] = I[m,k]/N[m]
            end
        end
    # end

    return a
end
