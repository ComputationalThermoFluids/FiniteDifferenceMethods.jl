(A::Type{<:SparseMatrixCSC})(this::ReshapedOperator) = convert(A, this)
(A::Type{SparseMatrixCSC})(this::ReshapedOperator) = convert(A, this)

_convert(::Type{SparseMatrixCSC}, op::Operator) =
    _convert(SparseMatrixCSC{eltype(op),Int}, op)
_convert(::Type{SparseMatrixCSC{T}}, op::Operator) where {T} =
    _convert(SparseMatrixCSC{T,Int}, op)

convert(A::Type{<:SparseMatrixCSC}, this::ReshapedOperator) =
    _convert(A, operator(this))

sparse(this::ReshapedOperator) = convert(SparseMatrixCSC, this)

convert(A::Type{<:SparseMatrixCSC}, this::Operator{T,1,1,2}) where {T} =
    _convert(A, this)
