"""

    spacing(n; a, b)

Computes spacing of a uniform mesh.


"""
spacing(n::Int; a=0, b=1) = (b-a) / (n-1)

"""

    collocated(n; a, b)

Creates uniform mesh of length `n` with end points matching `a` and `b`.


"""
collocated(n::Int; a=0, b=1) = range(a, b, step=spacing(n))
