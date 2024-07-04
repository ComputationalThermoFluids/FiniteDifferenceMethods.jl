"""

    OperatorSparsity


"""
abstract type OperatorSparsity end

struct HasStencil <: OperatorSparsity end
struct SparsityUnknown <: OperatorSparsity end

# default (operations are typically full)
OperatorSparsity(::Type{O}) where {O<:Operator} = SparsityUnknown()
#OperatorSparsity(::Type{O}) where {O<:∂} = HasStencil()

OperatorSparsity(::O) where {O<:Operator} = OperatorSparsity(O)

# interface

(op::Operator)(ind::CartInd, args::Varg{AArr}) =
    getcoef(OperatorSparsity(op), op, ind, args...)

getcoef(::SparsityUnknown, op::Operator, _...) = error("Specialization required.")
