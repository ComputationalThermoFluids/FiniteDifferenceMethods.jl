# something that can looks like a matrix but where rows and columns are multi-indexed
abstract type Operator{T,M,N,P} <: AArr{T,P} end

#
struct ReshapedOperator{T,A<:Operator{T}} <: AMat{T}
    op::A
end

# accessors

operator(this::ReshapedOperator) = this.op

# array interface

function size(this::ReshapedOperator)
    op = operator(this)
    length(rowinds(op)), length(colinds(op))
end

getindex(this::ReshapedOperator, i::Int, j::Int) = getindex(operator(this), i, j)

reshape(op::Operator, ::NColon{2}) = ReshapedOperator(op)
reshape(op::Operator, ::Vararg{Colon,2}) = ReshapedOperator(op)
