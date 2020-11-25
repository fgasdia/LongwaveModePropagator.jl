function romberg(x, y)
    # Adapted from NumericalIntegration.jl
    # This version always goes the maximum number of steps and uses MArrays for `rombaux`
    half = 1//2

    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    @assert ((length(x) - 1) & (length(x) - 2)) == 0 "Need length of vector to be 2^n + 1"
    length(x) < 2^20+1 || throw(ArgumentError("vector length must be less than 1048577"))

    maxsteps = convert(Int, log2(length(x)-1))

    rombaux = zeros(MMatrix{20,2,eltype(y),40})  # hardcoded to accomodate long x
    prevcol = 1
    currcol = 2
    @inbounds h = x[end] - x[1]
    @inbounds rombaux[1,1] = (y[1] + y[end])*h*half

    @inbounds for i in 1:(maxsteps-1)
        h *= half
        npoints = 1 << (i-1)
        jumpsize = div(length(x)-1, 2*npoints)
        c = zero(eltype(rombaux))
        for j in 1:npoints
            c += y[1 + (2*j-1)*jumpsize]
        end
        rombaux[1, currcol] = h*c + half*rombaux[1, prevcol]
        for j in 2:(i+1)
            n_k = exp2(j-1)^2  # ==  4^(j-1) but much faster!
            rombaux[j, currcol] = (n_k*rombaux[j-1, currcol] - rombaux[j-1, prevcol])/(n_k - 1)
        end

        prevcol, currcol = currcol, prevcol
    end

    @inbounds return rombaux[maxsteps, prevcol]
end
