using Plots; contourf, plot

struct WaveletCoherence <: AbstractArray{Complex{Float64}, 2}
    coh::Array{Complex{Float64}, 2}
    scales::Array{Float64, 1}
    wavelet::GenericContinuousWavelet
end

function coherence(wt::WaveletTransform, X::Array{Float64, 2})
    N,M = size(X)
    C = zeros(Complex{Float64}, N, length(wt.scales))
    P = zeros(Float64, N, length(wt.scales))
    for i in 1:M
        wX = transform(wt, X[:,i])
        C[:,:] += wX.wX
        P[:,:] += map(abs, wX.wX)
    end
    C ./= P

    WaveletCoherence(C, wt.scales, wt.wavelet)
end


Base.size(A::WaveletCoherence) = size(A.coh)
Base.getindex(A::WaveletCoherence, I::Vararg{Int, N}) where {N} = getindex(A.coh, I...);
Base.setindex!(A::WaveletCoherence, v, I::Vararg{Int, N}) where {N} = setindex!(A.coh, v, I...)

function Plots.contourf(A::WaveletCoherence; kw...)
    contourf(1:size(A.coh)[1], A.scales, transpose(map(abs, A.coh)); kw...);
    # draw valid range
end
