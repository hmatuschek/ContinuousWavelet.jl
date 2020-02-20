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

    threshold = undef;
    if nsurrogate > 0
        empC = Array{Float64, 1}[]
        # generate nsurrogate surrogate datasets and obtain coherence
        for i in 1:nsurrogate
            Cs = coherence(wt, surrogate(X); nsurrogate=0)
            append!(empC, reshape(Cs, sum(size(Cs))))
        end
        # get threshold corresponding to the given α-level
        sort!(empC)
        i = Int((1.0-α)*length(empC))
        thresholds = empC[(i<1) ? 1 : i]
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
    cplt = contourf(X, Y, Z; kw...);
    # TODO draw valid range
    # draw point-wise significance level if set
    if isdefined(A.α)
        contour!(X, Y, Y; levels=A.α, line_width=3, seriescolor=:reds)
    end
    return cplt
end
