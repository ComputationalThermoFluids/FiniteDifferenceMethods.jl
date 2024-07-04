"""

    Jacobian{N}(operator)

Jacobian of an operator with respect to a positional argument (`N`).

Provide following implementation for each differentiated operator:

```julia
(jac::∂)(ind::CartesianIndex, args::Vararg{AArr})

```

!!! note

    For linear dependencies, no need to specify: use `Undef`?


"""
struct Jacobian{N,O<:Operator} <: Operator
    op::O

    Jacobian{N}(op::O) where {N,O<:Operator} = new{N,O}(op)
end

const ∂ = Jacobian

Base.print_without_params(::Type{<:∂}) = false

#

operator(this::Jacobian) = this.op

#
#
#struct JacobianMatrix{T,P,M,N,O,R<:NTuple{M,AURange},C<:NTuple{N,AURange},A<:TupleN{AArr}} <: AMat{T}
#    op::O
#    rows::R
#    cols::C
#    args::A
#
#    function JacobianMatrix{P}(op::O, rows::R, cols::C, args::A) where {P,M,N,O<:Operator,R<:NTuple{M,AURange},C<:NTuple{N,AURange},A<:TupleN{AArr}}
#        T = Base.promote_eltype(args...)
#        new{T,P,M,N,O,R,C,A}(op, rows, cols, args)
#    end
#end
#
#const JacMat = JacobianMatrix
#
#operator(this::JacMat{T,P}) where {T,P} = ∂{P}(this.op)
#rowaxes(this::JacMat) = this.rows
#colaxes(this::JacMat) = this.cols
#arguments(this::JacMat) = this.args
#
##
#
#size(this::JacMat) = prod(length, rowaxes(this)), prod(length, colaxes(this))
#
#function getindex(this::JacMat{T}, i::Int, j::Int) where {T}
#    r = getindex(CartInds(rowaxes(this)), i)
#    c = getindex(CartInds(colaxes(this)), j)
#
#    convert(T, operator(this)(CartInd(r, c), arguments(this)...))
#end
#
## Linear algebra
#
#_lazyjac_to_mat(::Type{∂{1}}, op, rngs::NTuple{N,AURange}, x::ArrA{M}, args...) where {N,M} =
#    JacobianMatrix{1}(op, rngs[1:N-M], rngs[N-M+1:N], (x, args...))
#
#_lazy_reshape(this::∂{N}, rngs, args, ::Colon, ::Colon) where {N} =
#    _lazyjac_to_mat(∂{N}, operator(this), rngs, args...)
#
#reshape(this::Lazily, args::Vararg{Colon}) =
#    _lazy_reshape(operator(this), ranges(this), arguments(this), args...)

# Base.ReshapedArray{Float64, 1, Foo{Float64, 2}
_size(::∂{1}, rngs::NTup{N,AURange}, (x, _...)::Tup{ArrA{M},Varg{AArr}}, ::Varg{Colon,2}) where {N,M} =
    prod(length, rngs[1:N-M]),  prod(length, rngs[N-M+1:N])

# access by row

(op::∂)((i,)::Dims{1}, d::Tup{SInt{D}}, args...) where {D} = op(d, (i+D,), args...)
#(op::∂)(d::Tup{SInt{D}}, (j,)::Dims{1}, args...) where {D} = op((j-D,), d, args...)
