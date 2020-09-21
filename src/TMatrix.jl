"""
    TMatrix <: SMatrix{4, 4, ComplexF64, 16}

A custom `SMatrix` subtype that represents the `T` matrix from Clemmow and Heading
that is semi-sparse. The form of the matrix is:

```
┌ 11 12  0 14 ┐
| 0  0   1  0 |
| 31 32  0 34 |
└ 41 42  0 44 ┘
```

Implements custom (efficient) matrix-vector multiplication.
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

# TODO: Will this promote properly?
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
