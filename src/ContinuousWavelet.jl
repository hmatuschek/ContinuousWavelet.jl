module ContinuousWavelet
include("transform.jl")
include("surrogate.jl")
include("coherence.jl")
export CauchyWavelet, MorletWavelet, WaveletTransform, transform, surrogate
end  # module ContinuousWavelet
