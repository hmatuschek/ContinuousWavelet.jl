include("wavelet.jl")

using Plots; contourf, plot

@doc raw"""
Represents a wavelet transformed. Beside the wavelet coefficients, this object contains additional
information about the wavelet transformed including the scales at which the transform was performed
as well as the wavelet used. This object extends `AbstractArray{Complex{Float64}, 2}` thus allowing
for direct access to the wavelet coefficients. Each voice is stored as a column, thus the number of
rows matches the number of samples of the original time-series.

# Fields
- `wX`: The actual wavelet coefficients, each column represents a single voice.
- `scales`: The scales at which the wavelet transform was performed.
- `wavelet`: The wavelet of the transform.
"""
struct ContinuousWaveletTransformed <: AbstractArray{Complex{Float64}, 2}
    wX::Array{Complex{Float64}, 2}
    scales::Array{Float64, 1}
    wavelet::GenericContinuousWavelet
end

Base.size(A::ContinuousWaveletTransformed) = size(A.wX)
Base.getindex(A::ContinuousWaveletTransformed, I::Vararg{Int, N}) where {N} = getindex(A.wX, I...);
Base.setindex!(A::ContinuousWaveletTransformed, v, I::Vararg{Int, N}) where {N} = setindex!(A.wX, v, I...)

@doc raw"""
    Plots.contourf(A::ContinuousWavelet.ContinuousWaveletTransformed; kw...)

Plots the modulus (absolute value) of the wavelet transformed using a filled-contour plot.
Additional keyword arguments are passed to the default implementation of `Plots.contourf()`.
"""
function Plots.contourf(A::ContinuousWaveletTransformed; kw...)
    contourf(1:size(A.wX)[1], A.scales, transpose(map(abs, A.wX)); linewidth=0, kw...);
    # draw valid range
end
