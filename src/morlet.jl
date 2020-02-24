@doc raw"""
The Morlet wavelet.
```math
 g(t) = \sqrt{\frac{\delta}{2\pi}}\exp(2\pi\,i\,t-t^2\,\delta)\,,
```
where δ specifies the time-frequency resolution of the wavelet.
"""
struct MorletWavelet <: GenericContinuousWavelet
    dff::Float64
    norm::Float64

    function MorletWavelet(dff::Real)
        new(dff, dff/2/π);
    end
end

function eval_analysis(wav::MorletWavelet, t::Float64)
    exp(-t^2 * wav.dff/2 + 2π*t*1im) * wav.norm;
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
    # 99.9999% power at scale 1
    5/sqrt(wav.dff);
end;

function cutoff_freq(wav::MorletWavelet)
    # 99.9999% power at scale 1
    1 + 5*sqrt(wav.dff);
end;
