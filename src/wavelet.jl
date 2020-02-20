using SpecialFunctions: logabsgamma

abstract type GenericContinuousWavelet; end

include("morlet.jl")
include("cauchy.jl")
