using ContinuousWavelet: cutoff_time, eval_analysis

@testset "cutoff-time" begin
    @testset "Cauchy scale=$scale" for scale in [10, 20, 50, 100, 200, 500, 1000]
        wav = CauchyWavelet(10)
        t = Int32(ceil(cutoff_time(wav)*scale))
        P0 = sum(map(t -> abs(eval_analysis(wav, t/scale))^2, 0:t))
        P1 = sum(map(t -> abs(eval_analysis(wav, t/scale))^2, (t+1):(10*t)))
        @test P1/P0 < 1e-6
    end

    @testset "Morlet scale=$scale" for scale in [10, 20, 50, 100, 200, 500, 1000]
        wav = MorletWavelet(0.5)
        t = Int32(ceil(cutoff_time(wav)*scale))
        P0 = sum(map(t -> abs(eval_analysis(wav, t/scale))^2, 0:t))
        P1 = sum(map(t -> abs(eval_analysis(wav, t/scale))^2, (t+1):(10*t)))
        @test P1/P0 < 1e-6
    end
end
