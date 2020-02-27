@testset "synthesis" begin
    N  = 1024; t = 1:N; scale=100
    X  = sin.(2π*t/scale)
    wav = CauchyWavelet(10)
    wt = ContinuousWaveletTransform(wav, LinRange(2, 200, 64))
    wX = transform(wt, X)
    Y  = synthesis(wX)
    edge = Int32(ceil(scale*ContinuousWavelet.cutoff_time(wav)))
    @test maximum(abs.(X[edge:end-edge]-Y[edge:end-edge])) < 1e-2
end
