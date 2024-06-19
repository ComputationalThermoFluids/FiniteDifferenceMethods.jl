module FiniteDifferenceMethods

using LinearAlgebra,
      SparseArrays

import Base: convert,
             size,
             getindex

export spacing,
       collocated,
       staggered,
       laplacian,
       Gradient

include("aliases.jl")
include("utils.jl")
include("mesh.jl")
include("laplacian.jl")
include("gradient.jl")

end
