using SpecialFunctions: logabsgamma

abstract type GenericContinuousWavelet; end


struct MorletWavelet <: GenericContinuousWavelet
    dff::Float64
    norm::Float64

    function MorletWavelet(dff::Float64)
        new(dff, dff/2/π);
    end
end

function eval_analysis(wav::MorletWavelet, t::Float64)
    exp(-t^2 * wav.dff/2 + 2π*t*1im) * norm;
end

eval_synthesis(wav::MorletWavelet, t::Float64) = eval_analyis(wav, t);

function eval_repkern(wav::MorletWavelet, a::Float64, b::Float64)
    a2p1 = (a^2 + 1); am1 = (a - 1);
    d = wav.dff; d2 = d*d;
    c = d2/sqrt(2*π*d*a2p1);
    re = -(d2*b^2 + 4*π^2*am1^2)/(2*d*a2p1);
    im = 2*π*b*(1 + am1/a2p1)/a;
    c * exp(re + im*1im);
end

function cutoff_time(wav::MorletWavelet)
    # 99% power at scale 1
    3/sqrt(wav.dff);
end;

function cutoff_freq(wav::MorletWavelet)
    # 99.% power at scale 1
    1 + 3*sqrt(wav.dff);
end;


struct CauchyWavelet <: GenericContinuousWavelet
    α::Float64;
    norm::Float64;

    function CauchyWavelet(α::Float64)
        norm = exp(-2*log(2π) - logabsgamma(2α+1)[1]/2 + logabsgamma(α+1)[1]
          + (2α+1)*log(2)/2 + log(α));
        new(α, norm);
    end
end

function eval_analysis(wav::CauchyWavelet, t::Float64)
    (1 - 2im*π*t/wav.α)^(-1-wav.α);
end

eval_synthesis(wav::CauchyWavelet, t::Float64) = eval_analyis(wav, t);

function eval_repkern(wav::CauchyWavelet, a::Float64, b::Float64)
    c = wav.α*log(a) + logabsgamma(2*wav.α-1)[1] - (1+2*wav.α)*log(2π);
    exp(c) * (1+a-2im*π*b/wav.α)^(-1-2*wav.α);
end

function cutoff_time(wav::CauchyWavelet)
    eps = 1e-2;
    wav.α*sqrt(eps^(-2/(wav.α+1))-1)/(2π);
end;

function cutoff_freq(wav::CauchyWavelet)
    eps = 1e-2;
    1 + 1/( wav.α^2 * (eps^(-2 / (wav.α + 1)) - 1) / ((2π)^2) );
end;
