module ContinuousWavelet
include("transform.jl")
include("surrogate.jl")

export CauchyWavelet, MorletWavelet, WaveletTransform, transform, surrogate
end  # module ContinuousWavelet
