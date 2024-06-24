"""

`X`-contribution to volume-integrated divergence.

- `X`: component of the gradient (`:x`, `:y` or `:z`...),
- `V`: position in argument list (`:ω` or `:γ`).


"""
struct Divergence{X,V,M,T,A<:AArr{T,M},R<:TupleN{AURange},C<:TupleN{AURange},N} <: Operator{T,M,M,N}
    mom::NTuple{2,A}
    row::CartInds{M,R}
    col::CartInds{M,C}

    function Divergence{X,V}(mom::NTuple{2,A}, row::CartInds{M,R},
                             col::CartInds{M,C}) where {X,V,M,T,A<:AArr{T,M},R,C}
        new{X,V,M,T,A,R,C,2M}(mom, row, col)
    end
end

const Div = Divergence

# accessors

moments(this::Div) = this.mom
rowinds(this::Div) = this.row
colinds(this::Div) = this.col

# array interface

size(this::Div) = (size(rowinds(this))..., size(colinds(this))...)

function getindex(this::Div{X,V,M,T,A,R,C,N},
                  inds::Vararg{Int,N}) where {X,V,M,T,A,R,C,N}
    i, j = CartInd(inds[1:M]), CartInd(inds[M+1:N])
    getindex(this, i, j)
end

function getindex(this::Div, i::Int, j::Int)
    ci = getindex(eachindex(rowinds(this)), i)
    cj = getindex(eachindex(colinds(this)), j)
    getindex(this, ci, cj)
end

"""

!!! warning

    Assumes `i`*th* element of the gradient is collocated with ``x _ {i-1/2}``.


"""
function getindex(this::Div{:x,:ω,1,T}, i::CartInd{1}, j::CartInd{1}) where {T}
    A, _ = moments(this)

    r, = Tuple(getindex(rowinds(this), i))
    c, = Tuple(getindex(colinds(this), j))

    if isequal(c, r)
        -A[c] # qω[r]
    elseif isequal(c, r+1)
        +A[c] # qω[r+1]
    else
        zero(T)
    end
end

"""

!!! warning

    Assumes `i`*th* element of the gradient is collocated with ``x _ {i-1/2}``.


"""
function getindex(this::Div{:x,:γ,1,T}, i::CartInd{1}, j::CartInd{1}) where {T}
    A, B = moments(this)

    r, = Tuple(getindex(rowinds(this), i))
    c, = Tuple(getindex(colinds(this), j))

    if isequal(c, r)
        A[c] - B[r] # qγ[r]
    elseif isequal(c, r+1)
        B[r] - A[c] # qγ[r+1]
    else
        zero(T)
    end
end

getindex(this::Div, inds::CartInds) =
    getindex(this, inds.indices...)

function getindex(this::Div{X,V,M,T,A,R,C,N}, inds::Vararg{OrdinalRange,N}) where {X,V,M,T,A,R,C,N}
    row, col = CartInds(inds[1:M]), CartInds(inds[M+1:N])
    getindex(this, row, col)
end

getindex(this::Div{X,V,M}, row::CartInds{M}, col::CartInds{M}) where {X,V,M} =
    Div{X,V}(moments(this), getindex(rowinds(this), row), getindex(colinds(this), col))
