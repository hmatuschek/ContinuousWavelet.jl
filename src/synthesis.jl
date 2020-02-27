using DSP: conv

@doc raw"""
    synthesis(wt::WaveletTransformed)

Implements the continuous wavelet synthesis on some `WaveletTransformed` object. That is, it
reconstructs the signal, `wt` was created from. The optional keyword argument `hilbert` specifies
whether the complex Hilbert transformed of the original signal or the real signal is returned. By
default the real reconstruction is returned.
"""
function synthesis(wt::ContinuousWavelet.ContinuousWaveletTransformed; hilbert=false)
    N, M = size(wt)
    res = zeros(Complex{Float64}, N)
    # First convolution with smallest scale
    K = Int32(exp2(ceil(log2(2*cutoff_time(wt.wavelet)*wt.scales[1])))); off = K÷2
    kern = map(t -> eval_synthesis(wt.wavelet,
        (t-K/2)/wt.scales[1])/wt.scales[1]/wt.scales[1], 1:K)
    last = conv(kern, wt[:,1])[(off):(off+N-1)]

    current = similar(last)
    for i in 2:M
        # eval kernel for current scale
        K = Int32(exp2(ceil(log2(2*cutoff_time(wt.wavelet)*wt.scales[i])))); off = K÷2
        kern = map(t -> eval_synthesis(wt.wavelet,
            (t-K/2)/wt.scales[i])/wt.scales[i]/wt.scales[i], 1:K)
        # convolution of i-th voice with kernel at i-th scale
        current[:] = conv(kern, wt[:,i])[(off):(off+N-1)]
        # integrate by mid-point evaluation
        res[:] += ( (current+last) .* ((wt.scales[i]-wt.scales[i-1]))/2 )
        # swap current <-> last w/o copy
        tmp = last;
        last = current;
        current = tmp;
    end
    # If Hilbert transformed result is requested -> complex otherwise -> real
    (hilbert) ? res : 2*real.(res)
end
