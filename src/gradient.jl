"""

- `X`: component of the gradient (`:x`, `:y` or `:z`...),
- `V`: position in argument list (`:ω` or `:γ`).

!!! note

    Rename `row` `inner` and `col` `outer`?


"""
struct Gradient{X,V,M,T,A<:AArr{T,M},R<:TupleN{AURange},C<:TupleN{AURange},N} <: Operator{T,M,M,N}
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
                  ind::Vararg{Int,N}) where {X,V,M,T,A,R,C,N}
    i, j = CartInd(ind[1:M]), CartInd(ind[M+1:N])
    getindex(this, i, j)
end

function getindex(this::Grad, i::Int, j::Int)
    ci = getindex(eachindex(rowinds(this)), i)
    cj = getindex(eachindex(colinds(this)), j)
    getindex(this, ci, cj)
end

# diagonal coefficients

_main(op::Grad{X,:ω}, r, c) where {X} = ((_, B, W) = moments(op); +B[c] * pseudoinv(W[r]))
_lower(op::Grad{X,:ω}, r, c) where {X} = ((_, B, W) = moments(op); -B[c] * pseudoinv(W[r]))

_main(op::Grad{X,:γ}, r, c) where {X} = ((A, B, W) = moments(op); (A[r] - B[c]) * pseudoinv(W[r]))
_lower(op::Grad{X,:γ}, r, c) where {X} = ((A, B, W) = moments(op); (B[c] - A[r]) * pseudoinv(W[r]))

"""

!!! warning

    Assumes `i`*th* element of the gradient is collocated with ``x _ {i-1/2}``.


"""
function getindex(op::Grad{:x,V,1,T}, i::CartInd{1}, j::CartInd{1}) where {V,T}
    r = getindex(rowinds(op), i); rx, = Tuple(r)
    c = getindex(colinds(op), j); cx, = Tuple(c)

    if isequal(rx, cx)
        _main(op, r, c)
    elseif isequal(rx, cx+1)
        _lower(op, r, c)
    else
        zero(T)
    end
end

function getindex(op::Grad{:x,V,2,T}, i::CartInd{2}, j::CartInd{2}) where {V,T}
    r = getindex(rowinds(op), i); rx, ry = Tuple(r)
    c = getindex(colinds(op), j); cx, cy = Tuple(c)

    isequal(ry, cy) || return zero(T)

    if isequal(rx, cx)
        _main(op, r, c)
    elseif isequal(rx, cx+1)
        _lower(op, r, c)
    else
        zero(T)
    end
end

function getindex(op::Grad{:y,V,2,T}, i::CartInd{2}, j::CartInd{2}) where {V,T}
    r = getindex(rowinds(op), i); rx, ry = Tuple(r)
    c = getindex(colinds(op), j); cx, cy = Tuple(c)

    isequal(rx, cx) || return zero(T)

    if isequal(ry, cy)
        _main(op, r, c)
    elseif isequal(ry, cy+1)
        _lower(op, r, c)
    else
        zero(T)
    end
end

function getindex(op::Grad{:x,V,3,T}, i::CartInd{3}, j::CartInd{3}) where {V,T}
    r = getindex(rowinds(op), i); rx, ry, rz = Tuple(r)
    c = getindex(colinds(op), j); cx, cy, cz = Tuple(c)

    isequal(ry, cy) && isequal(rz, cz) || return zero(T)

    if isequal(rx, cx)
        _main(op, r, c)
    elseif isequal(rx, cx+1)
        _lower(op, r, c)
    else
        zero(T)
    end
end

function getindex(op::Grad{:y,V,3,T}, i::CartInd{3}, j::CartInd{3}) where {V,T}
    r = getindex(rowinds(op), i); rx, ry, rz = Tuple(r)
    c = getindex(colinds(op), j); cx, cy, cz = Tuple(c)

    isequal(rx, cx) && isequal(rz, cz) || return zero(T)

    if isequal(ry, cy)
        _main(op, r, c)
    elseif isequal(ry, cy+1)
        _lower(op, r, c)
    else
        zero(T)
    end
end

function getindex(op::Grad{:z,V,3,T}, i::CartInd{3}, j::CartInd{3}) where {V,T}
    r = getindex(rowinds(op), i); rx, ry, rz = Tuple(r)
    c = getindex(colinds(op), j); cx, cy, cz = Tuple(c)

    isequal(rx, cx) && isequal(ry, cy) || return zero(T)

    if isequal(rz, cz)
        _main(op, r, c)
    elseif isequal(rz, cz+1)
        _lower(op, r, c)
    else
        zero(T)
    end
end

getindex(this::Grad, inds::CartInds) =
    getindex(this, inds.indices...)

function getindex(this::Grad{X,V,M,T,A,R,C,N}, inds::Vararg{OrdinalRange,N}) where {X,V,M,T,A,R,C,N}
    row, col = CartInds(inds[1:M]), CartInds(inds[M+1:N])
    getindex(this, row, col)
end

getindex(this::Grad{X,V,M}, row::CartInds{M}, col::CartInds{M}) where {X,V,M} =
    Grad{X,V}(moments(this), getindex(rowinds(this), row), getindex(colinds(this), col))

# SparseArrays

"""

Leverage structure of backward difference form of the gradient operator.

!!! note

    Loops can easily be made thread-safe.


"""
function _convert(A::Type{<:SparseMatrixCSC{T,I}}, op::Grad{:x,V,1}) where {T,I,V}
    rs, = rowinds(op).indices
    cs, = colinds(op).indices

    # main diagonal
    dia = intersect(rs, cs)
    idia = unsafe_findin(rs, dia)
    jdia = unsafe_findin(cs, dia)

    # lower diagonal
    low = intersect(rs .- 1, cs)
    ilow = unsafe_findin(rs .- 1, low)
    jlow = unsafe_findin(cs, low)

    # for consistency with loop indexing with both contributions (*)
    if isempty(dia)
        idia = range(last(ilow)+1, length=0)
        jdia = range(last(jlow)+1, length=0)
    end

    m = length(rs)
    n = length(cs)
    p = length(dia) + length(low)

    colptr = Vector{I}(undef, n+1)
    rowval = Vector{I}(undef, p)
    nzval = Vector{T}(undef, p)

    j = firstindex(colptr)
    ptr = first(eachindex(rowval, nzval))

    # no contribution
    p = first(jlow) - firstindex(cs)
    colptr[j:j+p] .= ptr
    j += p

    bis = firstindex(ilow)

    # lower contribution only
    for c in cs[first(jlow):first(jdia)-1]
        r = c + 1

        i = ilow[bis]
        bis = nextind(ilow, bis)

        rowval[ptr] = i
        nzval[ptr] = _lower(op, r, c)
        ptr += 1

        j += 1
        colptr[j] = ptr
    end

    el = firstindex(idia)

    # both contributions (*)
    for c in cs[first(jdia):last(jlow)]
        # main
        r = c

        i = idia[el]
        el = nextind(idia, el)

        rowval[ptr] = i
        nzval[ptr] = _main(op, r, c)
        ptr += 1

        # lower
        r = c + 1

        i = ilow[bis]
        bis = nextind(ilow, bis)

        rowval[ptr] = i
        nzval[ptr] = _lower(op, r, c)
        ptr += 1

        j += 1
        colptr[j] = ptr
    end

    # main contribution only
    for c in cs[last(jlow)+1:last(jdia)]
        r = c

        i = idia[el]
        el = nextind(idia, el)

        rowval[ptr] = i
        nzval[ptr] = _main(op, r, c)
        ptr += 1

        j += 1
        colptr[j] = ptr
    end

    # no contribution
    p = lastindex(cs) - last(jdia)
    colptr[j+1:j+p] .= ptr
    j += p

    SparseMatrixCSC(m, n, colptr, rowval, nzval)
end
