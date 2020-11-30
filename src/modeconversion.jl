"""
this function takes just 1 mode conversion step
"""
function modeconversion(previous_wavefields::Wavefields{H},
                        wavefields::Wavefields{H},
                        adjoint_wavefields::Wavefields{H}) where H

    product = Vector{ComplexF64}(undef, numheights(wavefields))
    pproduct = similar(product)

    modes = eigenangles(wavefields)
    adjmodes = eigenangles(adjoint_wavefields)
    prevmodes = eigenangles(previous_wavefields)

    modes == adjmodes || @warn "Full mode conversion physics assumes adjoint
        wavefields are calculated for eigenangles of the original waveguide."

    nmodes = length(modes)
    nprevmodes = length(prevmodes)

    # `a` is total conversion
    a = Matrix{ComplexF64}(undef, nprevmodes, nmodes)
    for n in eachindex(modes)  # modes == adjmodes
        for m in eachindex(prevmodes)
            for i in eachindex(WAVEFIELD_HEIGHTS)
                @inbounds f = wavefields[i,n][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
                @inbounds fp = previous_wavefields[i,m][SVector(2,3,5,6)]  # Ey, Ez, Hy, Hz
                @inbounds g = SVector{4}(adjoint_wavefields[i,n][6],
                                         -adjoint_wavefields[i,n][5],
                                         -adjoint_wavefields[i,n][3],
                                         adjoint_wavefields[i,n][2])  # Hz, -Hy, -Ez, Ey

                gtranspose = transpose(g)
                product[i] = gtranspose*f
                pproduct[i] = gtranspose*fp
            end
            N = romberg(heights(wavefields), product)  # normalization
            I = romberg(heights(wavefields), pproduct)
            @inbounds a[m,n] = I/N
        end
    end

    return a
end
