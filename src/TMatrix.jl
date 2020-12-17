"""
    TMatrix <: SMatrix{4, 4, T}

A custom `SMatrix` subtype that represents the `T` matrix from [^Clemmow1954] that is
semi-sparse. The form of the matrix is:

```
┌                ┐
| T₁₁ T₁₂ 0  T₁₄ |
| 0   0   1  0   |
| T₃₁ T₃₂ 0  T₃₄ |
| T₄₁ T₄₂ 0  T₄₄ |
└                ┘
```

`TMatrix` implements efficient matrix-vector multiplication and other math based on its
special form.

See also: [`tmatrix`](@ref)

# References

[^Clemmow1954]: P. C. Clemmow and J. Heading, “Coupled forms of the differential equations
    governing radio propagation in the ionosphere,” Mathematical Proceedings of the
    Cambridge Philosophical Society, vol. 50, no. 2, pp. 319–333, Apr. 1954.
"""
struct TMatrix{T} <: StaticMatrix{4, 4, T}
    data::SVector{9,T}

    @inline function TMatrix(data::SVector{9,T}) where T
        new{T}(data)
    end
end

#==
# Constructors
==#

@generated function TMatrix(a::StaticMatrix{4, 4, T}) where {T}
    I = _tmatrix_indices()
    expr = Vector{Expr}(undef, 9)
    i = 0
    for index = 1:16
        if I[index] > 0
            expr[i += 1] = :(a[$index])
        end
    end
    quote
        @_inline_meta
        @inbounds return TMatrix(SVector{9,T}(tuple($(expr...))))
    end
end

@generated function TMatrix(a::NTuple{N,Any}) where {N}
    if N != 9
        throw(DimensionMismatch("Tuple length must be 9. Length of input was $N."))
    end
    expr = [:(a[$i]) for i = 1:9]
    quote
        @_inline_meta
        @inbounds return TMatrix(SVector{9}(tuple($(expr...))))
    end
end

#==
# Characteristics
==#

Base.size(::Type{TMatrix}) = (4, 4)
Base.size(::Type{TMatrix}, d::Int) = d > 2 ? 1 : 4

@inline function _tmatrix_indices()
    # Returns tuple containing TMatrix.data index for dense matrix index of the
    # tuple. For example, for linear index `5` of 4×4 matrix, return `4` which
    # is the index of TMatrix.data.
    # For "empty" values, return `0`.
    return tuple(1, 0, 2, 3, 4, 0, 5, 6, 0, 0, 0, 0, 7, 0, 8, 9)
end

@propagate_inbounds function Base.getindex(a::TMatrix{T}, i::Int) where T
    I = _tmatrix_indices()

    # Written out this is much faster than `i in (2, 6, 9, 11, 12, 14)`
    # Could probably be fewer lines with a generated function
    if I[i] > 0
        return a.data[I[i]]
    elseif i == 2
        return zero(T)
    elseif i ==  6
        return zero(T)
    elseif i == 9
        return zero(T)
    elseif i == 10
        return one(T)
    elseif i == 11
        return zero(T)
    elseif i == 12
        return zero(T)
    elseif i == 14
        return zero(T)
    else
        throw(BoundsError("attempt to access $(typeof(a)) at index $i"))
    end
end

@propagate_inbounds function getindex(a::TMatrix, inds::Int...)
    @boundscheck checkbounds(a, inds...)
    _getindex_scalar(a, inds...)
end

@generated function _getindex_scalar(a::TMatrix, inds::Int...)
    if length(inds) == 0
        return quote
            @_propagate_inbounds_meta
            a[1]
        end
    end

    # This function will only ever be called for length(inds) == 2
    stride = 1
    ind_expr = :()
    for i in 1:length(inds)
        if i == 1
            ind_expr = :(inds[1])
        else
            ind_expr = :($ind_expr + $stride * (inds[$i] - 1))
        end
        stride *= 4
    end
    return quote
        @_propagate_inbounds_meta
        a[$ind_expr]
    end
end

@generated function Base.Tuple(v::TMatrix)
    exprs = [:(v[$i]) for i = 1:16]
    quote
        @_inline_meta
        tuple($(exprs...))
    end
end

#==
# Math
==#

LinearAlgebra.ishermitian(v::TMatrix) = false
LinearAlgebra.issymmetric(v::TMatrix) = false

@inline Base.:(==)(a::TMatrix, b::TMatrix) = a.data == b.data

@inline Base.sum(v::TMatrix) = sum(v.data) + 1

@generated function Base.:*(a::Number, b::TMatrix{T}) where T
    I = _tmatrix_indices()
    exprs = [:(($I[$i] > 0) ? a*b.data[$I[$i]] : ($i == 10 ? convert(newtype,a) : zero(newtype))) for i = 1:16]
    quote
        @_inline_meta
        newtype = promote_type(T, typeof(a))
        SMatrix{4,4,newtype}(tuple($(exprs...)))
    end
end
@generated function Base.:*(a::TMatrix{T}, b::Number) where T
    I = _tmatrix_indices()
    exprs = [:(($I[$i] > 0) ? b*a.data[$I[$i]] : ($i == 10 ? convert(newtype,b) : zero(newtype))) for i = 1:16]
    quote
        @_inline_meta
        newtype = promote_type(T, typeof(b))
        SMatrix{4,4,newtype}(tuple($(exprs...)))
    end
end

@inline Base.:*(A::TMatrix, B::AbstractVector) = _mul(A, B)
@inline Base.:*(A::TMatrix, B::StaticVector) = _mul(Size(B), A, B)
@inline Base.:*(A::TMatrix, B::StaticMatrix) = _mul(Size(B), A, B)

@generated function _mul(a::TMatrix{Ta}, b::AbstractVector{Tb}) where {Ta,Tb}
    sa = Size(4,4)

    exprs = [:(a[$(LinearIndices(sa)[1,1])]*b[1]+a[$(LinearIndices(sa)[1,2])]*b[2]+a[$(LinearIndices(sa)[1,4])]*b[4]),
             :(b[3]),
             :(a[$(LinearIndices(sa)[3,1])]*b[1]+a[$(LinearIndices(sa)[3,2])]*b[2]+a[$(LinearIndices(sa)[3,4])]*b[4]),
             :(a[$(LinearIndices(sa)[4,1])]*b[1]+a[$(LinearIndices(sa)[4,2])]*b[2]+a[$(LinearIndices(sa)[4,4])]*b[4])]

    return quote
        @_inline_meta
        if length(b) != $sa[2]
            throw(DimensionMismatch("Tried to multiply arrays of size $sa and $(size(b))"))
        end
        T = Base.promote_op(LinearAlgebra.matprod, Ta, Tb)
        @inbounds return similar_type(b, T, Size($sa[1]))(tuple($(exprs...)))
    end
end

@generated function _mul(::Size{sb}, a::TMatrix{Ta}, b::StaticVector{<:Any, Tb}) where {sb, Ta, Tb}
    sa = Size(4,4)
    if sb[1] != sa[2]
        throw(DimensionMismatch("Tried to multiply arrays of size $sa and $sb"))
    end

    exprs = [:(a[$(LinearIndices(sa)[1,1])]*b[1]+a[$(LinearIndices(sa)[1,2])]*b[2]+a[$(LinearIndices(sa)[1,4])]*b[4]),
             :(b[3]),
             :(a[$(LinearIndices(sa)[3,1])]*b[1]+a[$(LinearIndices(sa)[3,2])]*b[2]+a[$(LinearIndices(sa)[3,4])]*b[4]),
             :(a[$(LinearIndices(sa)[4,1])]*b[1]+a[$(LinearIndices(sa)[4,2])]*b[2]+a[$(LinearIndices(sa)[4,4])]*b[4])]

    return quote
        @_inline_meta
        T = Base.promote_op(LinearAlgebra.matprod, Ta, Tb)
        @inbounds return similar_type(b, T, Size($sa[1]))(tuple($(exprs...)))
    end
end

@generated function _mul(Sb::Size{sb}, a::TMatrix{Ta}, b::StaticMatrix{<:Any,<:Any,Tb}) where {sb, Ta, Tb}
    sa = Size(4,4)
    if sb[1] != sa[2]
        throw(DimensionMismatch("Tried to multiply arrays of size $sa and $sb"))
    end

    S = Size(sa[1], sb[2])

    tmp_exprs = Vector{Vector{Expr}}(undef, sb[2])
    for i = 1:sb[2]
        tmp_exprs[i] = [:(a[$(LinearIndices(sa)[1,1])]*b[$(LinearIndices(sb)[1,i])]+
                          a[$(LinearIndices(sa)[1,2])]*b[$(LinearIndices(sb)[2,i])]+
                          a[$(LinearIndices(sa)[1,4])]*b[$(LinearIndices(sb)[4,i])]),
                        :(b[$(LinearIndices(sb)[3,i])]),
                        :(a[$(LinearIndices(sa)[3,1])]*b[$(LinearIndices(sb)[1,i])]+
                          a[$(LinearIndices(sa)[3,2])]*b[$(LinearIndices(sb)[2,i])]+
                          a[$(LinearIndices(sa)[3,4])]*b[$(LinearIndices(sb)[4,i])]),
                        :(a[$(LinearIndices(sa)[4,1])]*b[$(LinearIndices(sb)[1,i])]+
                          a[$(LinearIndices(sa)[4,2])]*b[$(LinearIndices(sb)[2,i])]+
                          a[$(LinearIndices(sa)[4,4])]*b[$(LinearIndices(sb)[4,i])])]
    end
    exprs = vcat(tmp_exprs...)

    return quote
        @_inline_meta
        T = Base.promote_op(LinearAlgebra.matprod, Ta, Tb)
        @inbounds return similar_type(a, T, $S)(tuple($(exprs...)))
    end
end

#==
Computing the T matrix entries
==#


@doc raw"""
    tmatrix(ea::EigenAngle, M)

Compute the matrix `T` as a `TMatrix` for the differential equations of a wave propagating
at angle `ea` in an ionosphere with susceptibility `M`.

Clemmow and Heading derived the `T` matrix from Maxwell's equations for an electromagnetic
wave in the anisotropic, collisional cold plasma of the ionosphere in a coordinate frame
where ``z`` is upward, propagation is directed obliquely in the ``x``-``z`` plane and
invariance is assumed in ``y``. For the four characteristic wave components
``e = (Ex, -Ey, Z₀Hx, Z₀Hy)ᵀ``, the differential equations are ``de/dz = -ikTe``.

See also: [`susceptibility`](@ref), [`dtmatrix`](@ref)

# References

[^Budden1955a]: K. G. Budden, “The numerical solution of differential equations governing
    reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227,
    no. 1171, pp. 516–537, Feb. 1955.

[^Clemmow1954]: P. C. Clemmow and J. Heading, “Coupled forms of the differential equations
    governing radio propagation in the ionosphere,” Mathematical Proceedings of the
    Cambridge Philosophical Society, vol. 50, no. 2, pp. 319–333, Apr. 1954.
"""
function tmatrix(ea::EigenAngle, M)
    S, C² = ea.sinθ, ea.cos²θ

    # Denominator of most of the entries of `T`
    den = inv(1 + M[3,3])

    M31den = M[3,1]*den
    M32den = M[3,2]*den

    T11 = -S*M31den
    T12 = S*M32den
    # T13 = 0
    T14 = (C² + M[3,3])*den
    # T21 = 0
    # T22 = 0
    # T23 = 1
    # T24 = 0
    T31 = M[2,3]*M31den - M[2,1]
    T32 = C² + M[2,2] - M[2,3]*M32den
    # T33 = 0
    T34 = S*M[2,3]*den
    T41 = 1 + M[1,1] - M[1,3]*M31den
    T42 = M[1,3]*M32den - M[1,2]
    # T43 = 0
    T44 = -S*M[1,3]*den

    # `TMatrix` is a special 4×4 matrix
    return TMatrix(T11, T31, T41,
                   T12, T32, T42,
                   T14, T34, T44)
end

"""
    dtmatrix(ea::EigenAngle, M)

Compute a dense `SMatrix` with the derivative of `T` with respect to `θ`.

See also: [`tmatrix`](@ref)
"""
function dtmatrix(ea::EigenAngle, M)
    S, C = ea.sinθ, ea.cosθ
    dC² = -2*S*C  # d/dθ (C²)

    den = inv(1 + M[3,3])

    dT11 = -C*M[3,1]*den
    dT12 = C*M[3,2]*den
    dT13 = 0
    dT14 = dC²*den
    dT21 = 0
    dT22 = 0
    dT23 = 0
    dT24 = 0
    dT31 = 0
    dT32 = dC²
    dT33 = 0
    dT34 = C*M[2,3]*den
    dT41 = 0
    dT42 = 0
    dT43 = 0
    dT44 = -C*M[1,3]*den

    return SMatrix{4,4}(dT11, dT21, dT31, dT41,
                        dT12, dT22, dT32, dT42,
                        dT13, dT23, dT33, dT43,
                        dT14, dT24, dT34, dT44)
end
