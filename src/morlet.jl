using SpecialFunctions: erfcinv

@doc raw"""
The Morlet wavelet.
```math
 g(t) = \sqrt{\frac{\delta}{2\pi}}\exp(2\pi\,i\,t-t^2\,\delta)\,,
```
where δ specifies the time-frequency resolution of the wavelet.
"""
struct MorletWavelet <: GenericContinuousWavelet
    δ::Float64
    norm::Float64
    σ::Float64

    @doc raw"""
        MorletWavelet(δ::Real; ϵ::Real=1e-3)

    Constructs a new Morlet wavelet, whith the given `dff` parameter specifying the time-frequency
    resolution. The optional keyword arguent ϵ specifies the cutoff at which the kernel evaluation
    gets truncated. It is defined as the fraction of total power loss of the mother wavelet. Smaller
    values of ϵ will increase the precision of the wavelet transform on the cost of longer kernels
    leading to slower convolutions.
    """
    function MorletWavelet(δ::Real; ϵ::Real=1e-3)
        new(δ, δ/2/π, sqrt(2)*erfcinv(2*ϵ));
    end
end

function eval_analysis(wav::MorletWavelet, t::Float64)
    exp(-t^2 * wav.δ/2 + 2π*t*1im) * wav.norm;
end

eval_synthesis(wav::MorletWavelet, t::Float64) = eval_analyis(wav, t);

function eval_repkern(wav::MorletWavelet, a::Float64, b::Float64)
    a2p1 = (a^2 + 1); am1 = (a - 1);
    d = wav.δ; d2 = d*d;
    c = d2/sqrt(2*π*d*a2p1);
    re = -(d2*b^2 + 4*π^2*am1^2)/(2*d*a2p1);
    im = 2*π*b*(1 + am1/a2p1)/a;
    c * exp(re + im*1im);
end

function cutoff_time(wav::MorletWavelet)
    # 99.9% power at scale 1
    wav.σ/sqrt(wav.δ);
end;

function cutoff_freq(wav::MorletWavelet)
    # 99.9% power at scale 1
    1 + wav.σ*sqrt(wav.δ);
end;
