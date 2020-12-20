"""
    modeconversion(previous_wavefields::Wavefields, wavefields::Wavefields,
        adjoint_wavefields::Wavefields; params=LMPParams())

Compute the mode conversion matrix `a` from the modes associated with `previous_wavefields`
to modes associated with `wavefields` and its `adjoint_wavefields`.

This is used in the approach known as full mode conversion [^Pappert1972b].

# References

[^Pappert1972b]: R. A. Pappert and R. R. Smith, “Orthogonality of VLF height gains in the
    earth ionosphere waveguide,” Radio Science, vol. 7, no. 2, pp. 275–278, 1972,
    doi: 10.1029/RS007i002p00275.
"""
function modeconversion(previous_wavefields::Wavefields{H}, wavefields::Wavefields{H},
                        adjoint_wavefields::Wavefields{H}; params=LMPParams()) where H

    @unpack wavefieldheights = params

    product = Vector{ComplexF64}(undef, numheights(wavefields))
    pproduct = similar(product)

    modes = eigenangles(wavefields)
    adjmodes = eigenangles(adjoint_wavefields)
    prevmodes = eigenangles(previous_wavefields)

    modes == adjmodes || @warn "Full mode conversion physics assumes adjoint wavefields are
        calculated for eigenangles of the original waveguide."

    nmodes = length(modes)
    nprevmodes = length(prevmodes)

    # `a` is total conversion
    a = Matrix{ComplexF64}(undef, nprevmodes, nmodes)
    for n in eachindex(modes)  # modes == adjmodes
        for m in eachindex(prevmodes)
            for i in eachindex(wavefieldheights)
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
            N, Nerr = romberg(heights(wavefields), product)  # normalization
            I, Ierr = romberg(heights(wavefields), pproduct)
            @inbounds a[m,n] = I/N
        end
    end

    return a
end
