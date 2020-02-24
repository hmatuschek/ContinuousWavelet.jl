@doc raw"""
The Cauchy or Paul wavelet.

In contrast to the `MorletWavelet`, the Cauchy or Paul wavelet is a proper wavelet with a
similar good localization in time and scale.

There are different definitions of the Cauchy wavelet around. Here one is implemented where
the center frequency is always 1 irrespective of the value of α:
```math
 g(t) = h(t) = \left(1-i\,\frac{2\pi\,t}{\alpha}\right)^{-(1+\alpha)}\,,
```
and its reproducing kernel
```math
 P_{g,h}(b, a) = \Gamma(2\alpha+1)\,a^{\alpha+1}\,\left(1+a-\frac{i\,b}{a}\right)^{-(2\alpha+1)}\,.
```
"""
struct CauchyWavelet <: GenericContinuousWavelet
    α::Float64;
    norm::Float64;

    function CauchyWavelet(α::Real)
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
    eps = 1e-3;
    wav.α*sqrt(eps^(-2/(wav.α+1))-1)/(2π);
end;

function cutoff_freq(wav::CauchyWavelet)
    eps = 1e-3;
    1 + 1/( wav.α^2 * (eps^(-2 / (wav.α + 1)) - 1) / ((2π)^2) );
end;
