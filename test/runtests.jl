using FiniteDifferenceMethods
using Test

@testset "FiniteDifferenceMethods.jl" begin
    n = 7

    A = laplacian((n,))
    x = collocated(n, start=0, stop=n+1)
    y = x .* (1 .- x)

    val = 2spacing(n+2) ^ 2
    @test all(isapprox(val), A * y)
end
