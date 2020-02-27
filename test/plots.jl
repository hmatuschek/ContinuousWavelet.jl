using Plots
using Random: randn;

@testset "plot wavelet transformed" begin
    X = randn(1024);
    wt = ContinuousWaveletTransform(CauchyWavelet(10), LinRange(2,200,64))
    wX = transform(wt, X)
    @test_nowarn contourf(wX)
end

@testset "plot coherence" begin
    N, M = 1024, 10;
    X = zeros(Float64, N, M)
    Y = zeros(Float64, N, M)
    for i in 1:M
        X[:,i] = map(t -> sin(2π*t/50), 1:N) + randn(N)
        Y[:,i] = map(t -> sin(2π*t/50), 1:N) + randn(N)
    end
    wt = ContinuousWaveletTransform(CauchyWavelet(10), LinRange(2,200,64))
    C = coherence(wt, X,Y)
    @test_nowarn contourf(C);
end
