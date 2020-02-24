module ContinuousWavelet
include("transform.jl")
include("surrogate.jl")
include("coherence.jl")
export CauchyWavelet, MorletWavelet, ContinuousWaveletTransform, transform, surrogate, surrogate!, coherence
end  # module ContinuousWavelet
