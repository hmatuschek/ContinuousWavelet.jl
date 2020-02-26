using Random: randn

@testset "coherence" begin
    N, M = 1024, 128
    t = 1:N
    X = zeros(Float64, N, M)
    for i in 1:M
        X[:,i] = sin.(2pi*t/50) + randn(N)
    end
    wt = ContinuousWaveletTransform(CauchyWavelet(10), [50])
    C = coherence(wt, X)
    @test sum(abs.(C)[:,1])/N ≈ 1.0 atol=1e-1
end
