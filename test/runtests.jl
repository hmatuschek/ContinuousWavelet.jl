using Test
using ContinuousWavelet

@testset "ContinuousWavelet.jl" begin
    include("surrogate.jl")
    include("cutofftime.jl")
    include("wtdelta.jl")
end
