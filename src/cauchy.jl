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
    ϵ::Float64

    @doc raw"""
        CauchyWavelet(α::Real; ϵ::Real=1e-2)

    Constructs a new Cauchy wavelet, whith the given α specifying the time-frequency resolution.
    The optional keyword arguent ϵ specifies the cutoff at which the kernel evaluation gets
    truncated. It is defined as the fraction of total power loss of the mother wavelet. Smaller
    values of ϵ will increase the precision of the wavelet transform on the cost of longer kernels
    leading to slower convolutions.
    """
    function CauchyWavelet(α::Real; ϵ::Real=1e-2)
        norm = exp(-2*log(2π) - logabsgamma(2α+1)[1]/2 + logabsgamma(α+1)[1]
          + (2α+1)*log(2)/2 + log(α));
        new(α, norm, ϵ);
    end
end

function eval_analysis(wav::CauchyWavelet, t::Float64)
    (1 - 2im*π*t/wav.α)^(-1-wav.α);
end

eval_synthesis(wav::CauchyWavelet, t::Float64) = eval_analysis(wav, t);

function eval_repkern(wav::CauchyWavelet, a::Float64, b::Float64)
    c = wav.α*log(a) + logabsgamma(2*wav.α-1)[1] - (1+2*wav.α)*log(2π);
    exp(c) * (1+a-2im*π*b/wav.α)^(-1-2*wav.α);
end

function cutoff_time(wav::CauchyWavelet)
    wav.α*sqrt(wav.ϵ^(-2/(wav.α+1))-1)/(2π);
end;

function cutoff_freq(wav::CauchyWavelet)
    1 + 1/( wav.α^2 * (wav.ϵ^(-2 / (wav.α + 1)) - 1) / ((2π)^2) );
end;
