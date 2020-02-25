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
    # Draw transformed as heatmap
    N,M = size(A)
    cn = contourf(1:N, A.scales, transpose(map(abs, A.wX));
         linewidth=0, linetype=:heatmap, xlims=(1,N), ylims=(minimum(A.scales),maximum(A.scales)),
         kw...);
    # draw valid range
    x = [1.0]; y = Float64[A.scales[1]]
    for i in 1:M
        tc = cutoff_time(A.wavelet)*A.scales[i]
        if (tc < N÷2)
            push!(y, A.scales[i])
            push!(x, tc)
        end
    end
    append!(x, reverse(N.-x))
    append!(y, reverse(y))
    plot!(x,y; linecolor=:black, linewidth=3, label="",
          fillrange=A.scales[end], fillalpha=0.5, fillcolor=:black)
    return cn
end
