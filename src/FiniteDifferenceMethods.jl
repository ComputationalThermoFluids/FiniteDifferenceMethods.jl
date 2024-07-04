module FiniteDifferenceMethods

using LinearAlgebra,
      SparseArrays,
      StaticArrays,
      Static

import Base: convert,
             size,
             axes,
             getindex,
             reshape

import SparseArrays: sparse

using StaticArrays: SUnitRange

export spacing,
       collocated,
       staggered,
       laplacian,
       ContributionStyle,
       Ω, Γ,
       CoordinateStyle,
       X, Y, Z,
       Operator,
       Lazily,
       operator,
       Jacobian,
       ∂,
#       JacobianMatrix,
       Gradient,
       LinearStencil,
       stencil

include("aliases.jl")
include("utils.jl")
include("mesh.jl")
include("laplacian.jl")
#include("arrays.jl")
#include("blocked.jl")
include("singletons.jl")
include("operator.jl")
include("jacobian.jl")
include("sparse.jl")
include("stencil.jl")
include("linear.jl")
include("gradient.jl")
#include("divergence.jl")
#include("reshape.jl")
#include("sparse.jl")
#include("extras.jl")

end
