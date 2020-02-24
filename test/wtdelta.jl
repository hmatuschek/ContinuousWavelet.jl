@testset "WT delta" begin
    N = 1024;
    scale = 100;

    X  = zeros(Float64, N); X[N÷2] = 1;
    Y  = map(t -> ContinuousWavelet.eval_analysis(CauchyWavelet(10), t/scale), LinRange(-N÷2+1,N÷2,N));

    wt = ContinuousWaveletTransform(CauchyWavelet(10), [scale]);
    wX = transform(wt, X);

    @test maximum(abs.(wX[:,1]-Y)) ≈ 0.0 atol=0.1
end
