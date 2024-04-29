"""

    Laplacian(dims, per)

Lazy `N`-dimensional Laplacian.


"""
struct Laplacian{N,T} <: AbstractMatrix{T}
    dims::Dims{N}
    per::NTuple{N,Bool}

    Laplacian{N,T}(dims, per) where {N,T} = new{N,T}(dims, per)
end

# constructor

Laplacian(dims::Dims{N}, ::Type{T}=Float64;
          per=ntuple(Returns(false), Val(N))) where {N,T<:Number} =
    Laplacian{N,T}(dims, per)

laplacian(dims::Dims{N}, ::Type{A}=SparseMatrixCSC{Float64,Int};
          per=ntuple(Returns(false), Val(N))) where {N,A<:AbstractMatrix} =
    convert(A, Laplacian(dims, eltype(A); per))


# accessor

shape(this::Laplacian) = this.dims
periodicity(this::Laplacian) = this.per

# array interface

size(this::Laplacian) = (prod(shape(this)), prod(shape(this)))

function getindex(this::Laplacian{N,T}, i::Int, j::Int) where {N,T}
    isequal(i, j) && return T(2N)

    dims, per = shape(this), periodicity(this)

    indices = CartesianIndices(dims)
    I, J = Tuple(indices[i]), Tuple(indices[j])

    n = sum(isequal.(I, J))
    isequal(n, N-1) || return zero(T)

    m = sum(neighbours.(dims, per, I, J))
    ifelse(isone(m), -one(T), zero(T))
end


"""

    neighbours(n, per, i, j)

Tests whether two points are closest neighbors.


"""
neighbours(n::Int, per::Bool, i::Int, j::Int) =
    isone(abs(i-j)) || (per && !isone(n) && isone(min(i, j)) && isequal(max(i, j), n))


"""

"""
function strang(n::Integer, isperiodic, ::Type{T}) where {T}
    if isperiodic
        if isone(n)
            spdiagm(0  => 2ones(T, n))
        elseif isequal(n, 2)
            spdiagm(-1  => -ones(T, n-1),
                     0  => 2ones(T, n),
                    +1  => -ones(T, n-1))
        else
            spdiagm(1-n => -ones(T, 1),
                    -1  => -ones(T, n-1),
                     0  => 2ones(T, n),
                    +1  => -ones(T, n-1),
                    n-1 => -ones(T, 1))
        end
    else
        if isone(n)
            spdiagm(0  => 2ones(T, n))
        else
            spdiagm(-1 => -ones(T, n-1),
                     0 => 2ones(T, n),
                    +1 => -ones(T, n-1))
        end
    end
end

_kron((Ix, Iy)::NTuple{2,AbstractMatrix}, (Jx, Jy)::NTuple{2,AbstractMatrix}) =
    kron(Iy, Jx) + kron(Jy, Ix)

_kron((Ix, Iy, Iz)::NTuple{3,AbstractMatrix}, (Jx, Jy, Jz)::NTuple{3,AbstractMatrix}) =
    kron(Iz, Iy, Jx) + kron(Iz, Jy, Ix) + kron(Jz, Iy, Ix)

function convert(::Type{SparseMatrixCSC{T,Int}}, this::Laplacian{1}) where {T}
    dims, per = shape(this), periodicity(this)
    strang(only(dims), only(per), T)
end

function convert(::Type{SparseMatrixCSC{T,Int}}, this::Laplacian{N}) where {T,N}
    dims, per = shape(this), periodicity(this)
    _kron(map(UniformScaling(one(T)), dims), strang.(dims, per, T))
end
