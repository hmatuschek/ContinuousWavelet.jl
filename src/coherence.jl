using Plots; contourf, plot

struct WaveletCoherence <: AbstractArray{Complex{Float64}, 2}
    coh::Array{Complex{Float64}, 2}
    scales::Array{Float64, 1}
    wavelet::GenericContinuousWavelet
    α::Float64
end

function coherence(wt::WaveletTransform, X::Array{Float64, 2}; nsurrogate=0, α=0.05)
    N,M = size(X)
    C = zeros(Complex{Float64}, N, length(wt.scales))
    P = zeros(Float64, N, length(wt.scales))
    for i in 1:M
        wX = transform(wt, X[:,i])
        C[:,:] += wX.wX
        P[:,:] += map(abs, wX.wX)
    end
    C ./= P

    threshold = NaN;
    if nsurrogate > 0
        empC = Array{Float64, 1}[]
        # generate nsurrogate surrogate datasets and obtain coherence
        for i in 1:nsurrogate
            Cs = coherence(wt, surrogate(X); nsurrogate=0).coh
            empC = vcat(empC, map(abs, reshape(Cs, (prod(size(Cs)),)) ) )
        end
        # get threshold corresponding to the given α-level
        sort!(empC)
        i = Int32(floor((1.0-α)*length(empC)))
        threshold = empC[(i<1) ? 1 : i]
    end
    WaveletCoherence(C, wt.scales, wt.wavelet, threshold)
end


Base.size(A::WaveletCoherence) = size(A.coh)
Base.getindex(A::WaveletCoherence, I::Vararg{Int, N}) where {N} = getindex(A.coh, I...);
Base.setindex!(A::WaveletCoherence, v, I::Vararg{Int, N}) where {N} = setindex!(A.coh, v, I...)

function Plots.contourf(A::WaveletCoherence; kw...)
    X = 1:size(A.coh)[1]
    Y = A.scales
    Z = transpose(map(abs, A.coh))
    cn = contourf(X, Y, Z; linewidth=0, kw...);
    # TODO draw valid range
    # draw point-wise significance level if set
    if isfinite(A.α)
        contour!(X, Y, Z; levels=[A.α], linewidth=3)
    end
    return cn
end
