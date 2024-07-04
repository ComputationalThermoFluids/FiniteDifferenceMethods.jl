"""

    ContributionStyle


"""
abstract type ContributionStyle end

struct BulkContribution <: ContributionStyle end

const Ω = BulkContribution

Base.print_without_params(::Type{Ω}) = false

struct BoundaryContribution <: ContributionStyle end

const Γ = BoundaryContribution

Base.print_without_params(::Type{Γ}) = false

"""

    CoordinateStyle


"""
abstract type CoordinateStyle end

struct CartesianCoordinate{N} <: CoordinateStyle end

const X = CartesianCoordinate{1}
const Y = CartesianCoordinate{2}
const Z = CartesianCoordinate{3}

"""

    ConcatenationStyle


"""
abstract type ConcatenationStyle end

struct VerticalConcatenation <: ConcatenationStyle end

# \vdots
const ⋮ = VerticalConcatenation

struct HorizontalConcatenation <: ConcatenationStyle end

# \cdots
const ⋯ = HorizontalConcatenation

"""

    StaggeringStyle


"""
abstract type StaggeringStyle end

struct ForwardStaggering <: StaggeringStyle end

# \blacktriangleright
const ▶ = ForwardStaggering

struct BackwardStaggering <: StaggeringStyle end

# \blacktriangleleft
const ◀ = BackwardStaggering
