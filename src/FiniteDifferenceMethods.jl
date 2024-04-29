module FiniteDifferenceMethods

using LinearAlgebra,
      SparseArrays

import Base: convert,
             size,
             getindex

export spacing,
       collocated,
       laplacian

include("mesh.jl")
include("laplacian.jl")

end
