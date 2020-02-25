module ContinuousWavelet
include("transform.jl")
include("surrogate.jl")
include("coherence.jl")
include("detrend.jl")
export CauchyWavelet, MorletWavelet, ContinuousWaveletTransform, transform, surrogate, surrogate!,
       coherence, detrend!, detrend
end  # module ContinuousWavelet
