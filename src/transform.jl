using DSP: conv

include("transformed.jl")

struct ContinuousWaveletTransform
    wavelet::GenericContinuousWavelet
    blocks::Array{Tuple{Int32,Int32,Int32}, 1}
    scales::Array{Float64, 1}
    kernels::Array{Complex{Float64}, 2}

    function ContinuousWaveletTransform(wav::GenericContinuousWavelet, scales::AbstractVector)
        blocks  = Tuple{Int32,Int32,Int32}[]
        lengths = Int32[]
        sscales = sort(scales)
        sizes   = Array{Int32}(exp2.(ceil.(log2.(cutoff_time(wav).*sscales))))
        kernels = zeros(Complex{Float64}, sizes[end-1], length(scales))

        # Start the first block for the first kernel
        nblock = 1
        push!(blocks, (1, 1, sizes[1]) )
        kernels[1:sizes[1],1] = map(t -> eval_analysis(wav, (t-sizes[1]/2)/sscales[1]), 1:sizes[1])
        # process the remaining kernels
        for i in 2:length(scales)
            kernels[1:sizes[i],i] = map(t -> eval_analysis(wav, (t-sizes[i]/2)/sscales[i]), 1:sizes[i])
            # if same size than last block:
            if blocks[end][2]==sizes[i]
                # -> add to block
                blocks[end] = (blocks[end][1], i, sizes[i]);
            else
                # otherwise add new block
                push!(blocks, (i,i,sizes[i]))
            end
        end

        new(wav, blocks, sscales, kernels)
    end
end

function transform(trans::ContinuousWaveletTransform, x::Array{Float64,1})
    N = length(x)
    M = length(trans.scales)
    wX = Array{Complex{Float64},2}(undef, N,M)
    for (a,b,K) in trans.blocks
        off = K÷2
        wX[:,(a:b)] = conv(trans.kernels[(1:K),(a:b)], x)[(off+1):(off+N)]
    end

    ContinuousWaveletTransformed(wX, trans.scales, trans.wavelet)
end
