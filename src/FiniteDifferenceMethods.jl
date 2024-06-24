module FiniteDifferenceMethods

using LinearAlgebra,
      SparseArrays,
      StaticArrays

import Base: convert,
             size,
             axes,
             getindex,
             reshape

import SparseArrays: sparse

export spacing,
       collocated,
       staggered,
       laplacian,
#       BlockedArray,
       Operator,
       Gradient,
       Divergence

include("aliases.jl")
include("utils.jl")
include("mesh.jl")
include("laplacian.jl")
#include("blocked.jl")
include("operator.jl")
include("gradient.jl")
include("divergence.jl")
#include("reshape.jl")
include("sparse.jl")

end
