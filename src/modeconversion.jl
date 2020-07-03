# infinity norm of the difference of the last two estimates
# basically `20` turns off warnings, but |real| values are typically >1e4
const ROMBERG_ACCURACY = 100.0

"""
this function takes just 1 mode conversion step

"""
function modeconversion(previous_wavefields::Wavefields{T},
                        wavefields::Wavefields{T}, adjwavefields::Wavefields{T}) where {T}

    @assert numheights(previous_wavefields) == numheights(wavefields)
    zs = heights(wavefields)

    # TODO: Assuming `length(zs)` is always the same, we can reuse `product` b/t calls
    product = Vector{ComplexF64}(undef, numheights(wavefields))

    modes = eigenangles(wavefields)
    adjmodes = eigenangles(adjwavefields)
    prevmodes = eigenangles(previous_wavefields)

    modes == adjmodes || @warn "Full mode conversion physics assumes adjoint
        wavefields are calculated for eigenangles of the original waveguide."

    nummodes = length(modes)
    numadjmodes = length(adjmodes)
    numprevmodes = length(prevmodes)

    # Calculate _inverse_ normalization terms, 1/N
    # the inverse is more efficient because we would otherwise repeatedly divide
    # by N later
    invN = Vector{ComplexF64}(undef, nummodes)
    for m in eachindex(modes)
        for i in eachindex(zs)
            @inbounds f = wavefields[m][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
            @inbounds g = SVector{4,ComplexF64}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
            product[i] = transpose(g)*f
        end

        # N = integrate(zs, product, RombergEven(ROMBERG_ACCURACY))
        N = romberg(zs, product)
        invN[m] = inv(N)
    end

    # TODO: because nummodes = numadjmodes, we don't need to have two sets of loops
    # we can just have a scalar N, but will need to flip dimensions of `a` to
    # write down column first. Also, we would then need two product vectors, which
    # are relatively large

    # `a` is total conversion
    a = Matrix{ComplexF64}(undef, numadjmodes, numprevmodes)
    for k in eachindex(prevmodes)
        for m in eachindex(adjmodes)
            for i in eachindex(zs)
                @inbounds f = previous_wavefields[k][i][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
                @inbounds g = SVector{4,ComplexF64}(adjwavefields[m][i][6], -adjwavefields[m][i][5], -adjwavefields[m][i][3], adjwavefields[m][i][2])  # Hz, -Hy, -Ez, Ey
                product[i] = transpose(g)*f
            end
            # I = integrate(zs, product, RombergEven(ROMBERG_ACCURACY))
            I = romberg(zs, product)
            @inbounds a[m,k] = I*invN[m]
        end
    end

    return a
end
