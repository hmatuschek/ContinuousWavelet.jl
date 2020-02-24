using ContinuousWavelet
using Test
using Random


@testset "surrogate" begin
    N, M = 1024, 2
    X = zeros(Float64, N, M)
    for i in 1:M
        X[:,i] = map(n -> sin(2π*n/200), 1:N)
    end

    # Test 1D surrogate norm conservation
    Y = surrogate(X[:,1])
    normX = sum(map(t -> abs(t)^2, X[:,1]))
    normY = sum(map(t -> abs(t)^2, Y))
    @test normX ≈ normY rtol=1/N

    # Test 1D in-place surrogate norm conservation
    Y = copy(X[:,1]); surrogate!(Y)
    normX = sum(map(t -> abs(t)^2, X[:,1]))
    normY = sum(map(t -> abs(t)^2, Y))
    @test normX ≈ normY rtol=1/N

    # Test 2D surrogate norm conservation
    Y = surrogate(X)
    normX = sum(map(t -> abs(t)^2, X))
    normY = sum(map(t -> abs(t)^2, Y))
    @test normX ≈ normY rtol=1/N

    # Test 2D in-place surrogate norm conservation
    Y = copy(X); surrogate!(Y)
    normX = sum(map(t -> abs(t)^2, X))
    normY = sum(map(t -> abs(t)^2, Y))
    @test normX ≈ normY rtol=1/N
end
