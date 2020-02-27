using Test
using ContinuousWavelet

@testset "ContinuousWavelet.jl" begin
    include("detrend.jl")
    include("surrogate.jl")
    include("cutofftime.jl")
    include("wtdelta.jl")
    include("coherence.jl")
    include("synthesis.jl")
    include("plots.jl")
end
