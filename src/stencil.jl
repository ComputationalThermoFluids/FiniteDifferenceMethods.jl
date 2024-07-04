"""

    abstract type Stencil end

!!! note

    Rename `SparsityPattern`?


"""
abstract type Stencil end

struct LinearStencil{X,R<:SUnitRange} <: Stencil
    rng::R

    LinearStencil{X}(rng::R) where {X,R<:SUnitRange} =
        new{X,R}(rng)
end

# interface

getcoef(::HasStencil, op, ind, args...) =
    getstencil(stencil(op), op, ind, args...)

using Base.Cartesian

@generated function getstencil(::LinearStencil{X,SURange{S,L}}, op::∂{1},
                               ind::CartInd{2}, x::AVec, args::AArr...) where {S,L}
    quote
        i, j = Tuple(ind)
        k = j - i

        @nif($(L+1), d -> isequal(d+$S-1, k),
#                     d -> op((i,), (SInt{k}(),), x, args...),
                     d -> op((SInt{k}(),), (j,), x, args...),
                     d -> zero(Base.promote_eltype(x, args...)))
    end
end
