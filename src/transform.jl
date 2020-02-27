using DSP: conv

include("transformed.jl")

@doc raw"""
This struct represents a planned continuous wavelet transform. It contains the voices on which the
transform is performed as well the pre-computed analysis wavelets.

# Fields
- `wavelet`: Holds a reference to the wavelet to use for the transformation.
- `blocks`: The grouping of the voices of same kernel size for faster block-convolution.
- `scales`: The scales on which the wavelet transform will be performed.
- `kernels`: Matrix of evaluated kernels for the convolution.
"""
struct ContinuousWaveletTransform
    wavelet::GenericContinuousWavelet
    blocks::Array{Tuple{Int32,Int32,Int32}, 1}
    scales::Array{Float64, 1}
    kernels::Array{Complex{Float64}, 2}

    @doc raw"""
        ContinuousWaveletTransform(wav::GenericContinuousWavelet, scales::AbstractVector)

    Constructs a new wavelet transform for the given wavelet (`wav`) and `scales`. This constructor
    pre-computes the kernels for the given scales. It is therefore efficient to perform several
    transformations using `transform()` with the same `ContinuousWaveletTransform` object.
    """
    function ContinuousWaveletTransform(wav::GenericContinuousWavelet, scales::AbstractVector)
        blocks  = Tuple{Int32,Int32,Int32}[]
        lengths = Int32[]
        sscales = sort(scales)
        sizes   = Array{Int32}(exp2.(ceil.(log2.(2*cutoff_time(wav).*sscales))))
        kernels = zeros(Complex{Float64}, sizes[end], length(scales))

        # Start the first block for the first kernel
        push!(blocks, (1, 1, sizes[1]) )
        kernels[1:sizes[1],1] = map(t -> eval_analysis(wav, (t-sizes[1]/2)/sscales[1])/sscales[1], 1:sizes[1])
        # process the remaining kernels
        for i in 2:length(scales)
            kernels[1:sizes[i],i] = map(t -> eval_analysis(wav, (t-sizes[i]/2)/sscales[i])/sscales[i], 1:sizes[i])
            # if same size than last block:
            if blocks[end][3]==sizes[i]
                # -> add to block
                blocks[end] = (blocks[end][1], i, sizes[i]);
            else
                # otherwise add new block
                push!(blocks, (i,i,sizes[i]))
            end
        end
        nBlocks = length(blocks)
        minBlockSize = minimum(sizes)
        maxBlockSize = maximum(sizes)
        new(wav, blocks, sscales, kernels)
    end
end

@doc raw"""
    transform(trans::ContinuousWaveletTransform, x::AbstractArray{Float64,1})

Performes the planned continuous wavelet transfrom given by `trans` of the timeseries `x`. The
result is a `ContinuousWavelet.ContinuousWaveletTransformed` object.
"""
function transform(trans::ContinuousWaveletTransform, x::AbstractArray{Float64,1})
    N = length(x)
    M = length(trans.scales)
    wX = Array{Complex{Float64},2}(undef, N,M)
    for (a,b,K) in trans.blocks
        off = K÷2
        wX[:,(a:b)] = conv(trans.kernels[(1:K),(a:b)], x)[(off):(off+N-1)]
    end

    ContinuousWaveletTransformed(wX, trans.scales, trans.wavelet)
end
