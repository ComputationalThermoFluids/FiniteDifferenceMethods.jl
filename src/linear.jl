function sparse(this::Lazily)
    op = operator(this)
    _sparse(OperatorSparsity(op), this)
end

_sparse(::SparsityUnknown, this::Lazily{T,2}) where {T} =
    convert(SparseMatrixCSC, this)

function _sparse(::HasStencil, this::Lazily{T,2}) where {T}
    op = operator(this)
    rngs = ranges(this)
    args = arguments(this)

    __sparse(stencil(op), op, rngs, args)
end

function __sparse(::LinearStencil{X,SURange{-1,2}}, op::∂{1},
                  rngs::NTup{2,AURange}, args::Tup{ArrA{1},Varg{AArr}})
    T = Base.promote_eltype(args...)

    rs, cs = rngs

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

    colptr = Vector{Int}(undef, n+1)
    rowval = Vector{Int}(undef, p)
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
#        nzval[ptr] = _lower(op, r, c)
#        nzval[ptr] = op((r,), (SInt{-1}(),), args...)
        nzval[ptr] = op((SInt{-1}(),), (c,), args...)
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
#        nzval[ptr] = _main(op, r, c)
#        nzval[ptr] = op((r,), (SInt{0}(),), args...)
        nzval[ptr] = op((SInt{0}(),), (c,), args...)
        ptr += 1

        # lower
        r = c + 1

        i = ilow[bis]
        bis = nextind(ilow, bis)

        rowval[ptr] = i
#        nzval[ptr] = _lower(op, r, c)
#        nzval[ptr] = op((r,), (SInt{-1}(),), args...)
        nzval[ptr] = op((SInt{-1}(),), (c,), args...)
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
#        nzval[ptr] = _main(op, r, c)
#        nzval[ptr] = op((r,), (SInt{0}(),), args...)
        nzval[ptr] = op((SInt{0}(),), (c,), args...)
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
