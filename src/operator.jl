"""

An operator is a function that takes `P` `AbstractArray{T,N}`s and returns another `AbstractArray{T,M}`.


"""
abstract type Operator <: Function end


"""

Lazy representation of an operation (*i.e.* the application of an operator to a set of arguments).

!!! note

    `LazyOperation` indexing is **one-based**.

"""
struct LazyOperation{T,N,O<:Operator,R<:NTuple{N,AURange},A<:TupleN{AArr}} <: AArr{T,N}
    op::O
    rngs::R
    args::A

    function LazyOperation(op::O, rngs::R, args::A) where {N,O<:Operator,R<:NTuple{N,AURange},A<:TupleN{AArr}}
        T = Base.promote_eltype(args...)
        new{T,N,O,R,A}(op, rngs, args)
    end
end

const Lazily = LazyOperation

Base.print_without_params(::Type{<:Lazily}) = false

# accessors

operator(this::Lazily) = this.op
ranges(this::Lazily) = this.rngs
arguments(this::Lazily) = this.args

# array interface

size(op::Lazily) = length.(ranges(op))

function getindex(this::Lazily{T,N}, I::Vararg{Int,N}) where {T,N}
    ind = CartesianIndex(getindex.(ranges(this), I))
    convert(T, operator(this)(ind, arguments(this)...))
end

# out-of-place defaults to lazy representation
(op::Operator)(rngs::TupleN{AURange}, args::Vararg{AArr}) =
    Lazily(op, rngs, args)

# in-place requires specialization
 (op::Operator)(y::AArr, rngs::TupleN{AURange}, args::Vararg{AArr}) =
    (warn(lazy"In-place $(typeof(op)) requires specialization."); y)

# reshape

reshape(this::Lazily, args::Colon...) =
    reshape(this, _size(operator(this), ranges(this), arguments(this), args...))

reshape(this::Lazily{T,2}, ::Varg{Colon,2}) where {T} = this
