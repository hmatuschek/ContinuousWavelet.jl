include("wavelet.jl")

using Plots; contourf, plot

struct WaveletTransformed <: AbstractArray{Complex{Float64}, 2}
    wX::Array{Complex{Float64}, 2}
    scales::Array{Float64, 1}
    wavelet::GenericContinuousWavelet
end

Base.size(A::WaveletTransformed) = size(A.wX)
Base.getindex(A::WaveletTransformed, I::Vararg{Int, N}) where {N} = getindex(A.wX, I...);
Base.setindex!(A::WaveletTransformed, v, I::Vararg{Int, N}) where {N} = setindex!(A.wX, v, I...)

function Plots.contourf(A::WaveletTransformed; kw...)
    contourf(1:size(A.wX)[1], A.scales, transpose(map(abs, A.wX)); linewidth=0, kw...);
    # draw valid range
end
