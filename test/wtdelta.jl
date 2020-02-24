using ContinuousWavelet
using Test

@testset "WT delta" begin
    N = 128*1024;
    scale = 200;

    X  = zeros(Float64, N);
    X[N÷2] = 1;
    wt = ContinuousWaveletTransform(CauchyWavelet(10.0), [scale]);
    wX = transform(wt, X);

    print(maximum(abs.(wX[:,1] -
          map(t -> ContinuousWavelet.eval_analysis(CauchyWavelet(10.0), t/scale),
              LinRange(-N÷2+1,N÷2,N)))))
    @test wX[:,1] ≈ map(t -> ContinuousWavelet.eval_analysis(CauchyWavelet(10.0), t/scale),
                        LinRange(-N÷2+1,N÷2,N)) atol=1e-5
end
