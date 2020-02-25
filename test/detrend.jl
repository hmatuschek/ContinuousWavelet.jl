using Random: randn

@testset "detrend()" begin
    N = 1000
    X = randn(N, 2)
    Y = copy(X)
    for i in 1:N
        Y[i,:] = X[i,:] .+ 1 .+ (2*i)
    end

    @test X ≈ detrend(Y) rtol=3/sqrt(N);

    detrend!(Y)
    @test X ≈ Y rtol=3/sqrt(N)
end
