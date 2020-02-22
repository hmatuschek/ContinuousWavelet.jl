using SpecialFunctions: logabsgamma

@doc raw"""
Represents the base-class of all continuous wavelets.
"""
abstract type GenericContinuousWavelet; end

@doc raw"""
    eval_analysis(wav::ContinuousWavelet.GenericContinuousWavelet, t::Float64)

Evaluates the *mother wavelet* at timepoint `t`.
"""
function eval_analysis end

@doc raw"""
    eval_synthesis(wav::ContinuousWavelet.GenericContinuousWavelet, t::Float64)

Evaluates the *mother synthesis wavelet* at timepoint `t`.
"""
function eval_synthesis end

@doc raw"""
    eval_repkern(wav::ContinuousWavelet.GenericContinuousWavelet, a::Float64, b::Float64)

Evaluates the reproducing kernel at scale `a` and time-point `b`.
"""
function eval_repkern end

@doc raw"""
    cutoff_time(wav::ContinuousWavelet.GenericContinuousWavelet)

Returns the cut-off time-point for the given wavelet at scale 1.
"""
function cutoff_time end

@doc raw"""
    cutoff_freq(wav::ContinuousWavelet.GenericContinuousWavelet)

Returns the cut-off frequency for the given wavelet at scale 1.
"""
function cutoff_freq end


include("morlet.jl")
include("cauchy.jl")
