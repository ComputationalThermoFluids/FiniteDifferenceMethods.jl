"""

- `X`: component of the gradient (`:x`, `:y` or `:z`...),
- `V`: position in argument list (`:ω` or `:γ`).

"""
struct Gradient{X,V,M,T,A<:AArr{T,M},R,C,N} <: AArr{T,N}
    mom::NTuple{3,A}
    row::CartInds{M,R}
    col::CartInds{M,C}

    function Gradient{X,V}(mom::NTuple{3,A}, row::CartInds{M,R},
                           col::CartInds{M,C}) where {X,V,M,T,A<:AArr{T,M},R,C}
        new{X,V,M,T,A,R,C,2M}(mom, row, col)
    end
end

const Grad = Gradient

# accessors

moments(this::Grad) = this.mom
rowinds(this::Grad) = this.row
colinds(this::Grad) = this.col

# array interface

size(this::Grad) = (size(rowinds(this))..., size(colinds(this))...)

function getindex(this::Grad{X,V,M,T,A,R,C,N},
                  inds::Vararg{Int,N}) where {X,V,M,T,A,R,C,N}
    i, j = CartInd(inds[1:M]), CartInd(inds[M+1:N])
    getindex(this, i, j)
end

"""

!!! warning

    Assumes `i`*th* element of the gradient is collocated with ``x _ {i-1/2}``.


"""
function getindex(this::Grad{:x,:ω,1,T}, i::CartInd{1}, j::CartInd{1}) where {T}
    A, B, W = moments(this)

    r, = Tuple(getindex(rowinds(this), i))
    c, = Tuple(getindex(colinds(this), j))

    if isequal(c, r)
        +B[c] * pseudoinv(W[r])
    elseif isequal(c, r-1)
        -B[c] * pseudoinv(W[r])
    else
        zero(T)
    end
end

"""

!!! warning

    Assumes `i`*th* element of the gradient is collocated with ``x _ {i-1/2}``.


"""
function getindex(this::Grad{:x,:γ,1,T}, i::CartInd{1}, j::CartInd{1}) where {T}
    A, B, W = moments(this)

    r, = Tuple(getindex(rowinds(this), i))
    c, = Tuple(getindex(colinds(this), j))

    if isequal(c, r)
        (A[r] - B[c]) * pseudoinv(W[r])
    elseif isequal(c, r-1)
        (B[c] - A[r]) * pseudoinv(W[r])
    else
        zero(T)
    end
end

getindex(this::Grad, inds) =
    getindex(this, inds.indices...)

function getindex(this::Grad{X,V,M,T,A,R,C,N}, inds::Vararg{OrdinalRange{Int},N}) where {X,V,M,T,A,R,C,N}
    row, col = CartInds(inds[1:M]), CartInds(inds[M+1:N])
    getindex(this, row, col)
end

getindex(this::Grad{X,V,M}, row::CartInds{M}, col::CartInds{M}) where {X,V,M} =
    Grad{X,V}(moments(this), getindex(rowinds(this), row), getindex(colinds(this), col))
