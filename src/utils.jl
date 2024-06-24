ncolon(v::Val) = ntuple(Returns(:), v)

pseudoinv(x) = ifelse(iszero(x), x, inv(x))

"""

    unsafe_findin(a, b)

Return the range of indices in `a` where values are also in `b` assuming `b` is a subset of `a`.


"""
function unsafe_findin(a::AURange, b::AURange)
    @assert issubset(b, a)

    if isempty(b)
        start = firstindex(a)
        stop = firstindex(a)-1
    else
        start = searchsortedfirst(a, first(b))
        stop = searchsortedlast(a, last(b))
    end

    range(start, stop)
end
