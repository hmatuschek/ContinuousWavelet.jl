using Plots; contourf, plot

@doc raw"""
Represents a continuous wavelet coherence analysis.

# Fields:
- `coh`: The actual coherence as a complex value. The layout of this array follows the layout of
  `ContinuousWaveletTransformed`. That is, each voice is stored in a column, starting with the
  smallest scale.
- `scales`: The vector of scales at which the wavelet transform is performed.
- `wavelet`: The wavelet being used.
- `α`: The point-wise empirical signigicance level of the coherence analysis obtained by means of
  surrogate data. Holds `NaN` if no skipped.
"""
struct ContinuousWaveletCoherence <: AbstractArray{Complex{Float64}, 2}
    coh::Array{Complex{Float64}, 2}
    scales::Array{Float64, 1}
    wavelet::GenericContinuousWavelet
    α::Float64
end

@doc raw"""
    coherence(wt::ContinuousWaveletTransform,
              X::AbstractArray{Float64, 2};
              nsurrogate=0, α=0.05)

Implements the *self-coherence* of a set of time-series. The *self-coherence* is defined as
```math
 C(X)(b,a) = \frac{\left\vert E\left[\mathcal{W}[X](b,a)\right]\right\vert^2}{E\left[\left\vert\mathcal{W}[X](b,a)\right\vert^2\right]} =
  \frac{\left\vert\mathcal{W}[E\left[X\right]](b,a)\right\vert^2}{E\left[\left\vert\mathcal{W}[X](b,a)\right\vert^2\right]}
```
This measure can be interpreted as how coherent a process is from realization to realization.
Typical applications are the analyses of EEG event-related potentials and similar time-series of
repeated experiments.

# Arguments
- `wt`: Specifies the wavelet transform to use.
- `X`: Specifies the time-series to analyze. Each single time-series (trial) is stored as a column in `X`.
- 'nsurrogate': A keyword argument specifying the number of surrogate samples to generate from the
  given time-series for the estimation of the point-wise significant test. By default it is 0 and
  no surrogate samples are generated.
- `α`: A keyword arguement specifying the significance level for the point-wise significance test.
  By default set to `α=0.05`.
"""
function coherence(wt::ContinuousWaveletTransform, X::AbstractArray{Float64, 2}; nsurrogate=0, α=0.05)
    N,M = size(X)
    C = zeros(Complex{Float64}, N, length(wt.scales))
    P = zeros(Float64, N, length(wt.scales))
    for i in 1:M
        wX = transform(wt, X[:,i])
        C[:,:] += wX.wX
        P[:,:] += map(z->z^2, map(abs, wX.wX))
    end
    C = map(z->z^2, map(abs, C)) ./ M;
    C ./= P


    threshold = NaN;
    if nsurrogate > 0
        empC = Array{Float64, 1}[]
        # generate nsurrogate surrogate datasets and obtain coherence
        for i in 1:nsurrogate
            Cs = coherence(wt, surrogate(X); nsurrogate=0).coh
            empC = vcat(empC, map(abs, reshape(Cs, (prod(size(Cs)),)) ) )
        end
        # get threshold corresponding to the given α-level
        sort!(empC)
        i = Int32(floor((1.0-α)*length(empC)))
        threshold = empC[(i<1) ? 1 : i]
    end
    ContinuousWaveletCoherence(C, wt.scales, wt.wavelet, threshold)
end

@doc raw"""
    coherence(wt::ContinuousWaveletTransform,
              X::AbstractArray{Float64, 2}, Y::AbstractArray{Float64, 2};
              nsurrogate=0, α=0.05)

Canonical wavelet coherence as defined by
```math
 C(X)(b,a) = \frac{E\left[\mathcal{W}[X](b,a)\,\mathcal{W}[Y]^{*}(b,a)\right]}{
  E\left[\left\vert\mathcal{W}[X](b,a)\right\vert\,\left\vert\mathcal{W}[Y](b,a)\right\vert\right]}
```
That is the normalized expectation value of the product of the two wavelet transformed, where
`z^*` is the complex conjugate.

In contrast to the *self-coherence*, this analysis does estimate the coherence of a process to
itself but the coherence between two observed processes. It is therefore necessary that the two
sets of time-series `X` and `Y` are of the same shape. That is, each trial of `X` and `Y` must
consists of the same number of samples (rows) and both sets must contain the same number of trials
(columns).

# Arguments
- `wt`: Specifies the wavelet transform to use.
- `X`: Specifies one time-series to analyze. Each single time-series (trial) is stored as a column
    in `X`. Must be of same shape as `Y`.
- `Y`: Specifies one time-series to analyze. Each single time-series (trial) is stored as a column
    in `Y`. Must be of same shape as `X`.
- 'nsurrogate': A keyword argument specifying the number of surrogate samples to generate from the
    given time-series for the estimation of the point-wise significant test. By default it is 0
    and no surrogate samples are generated.
- `α`: A keyword arguement specifying the significance level for the point-wise significance test.
    By default set to `α=0.05`.
"""
function coherence(wt::ContinuousWaveletTransform, X::AbstractArray{Float64, 2}, Y::AbstractArray{Float64, 2}; nsurrogate=0, α=0.05)
    (size(X) != size(Y)) && error("Arrays X & Y must be of same size!");
    N,M = size(X)
    C = zeros(Complex{Float64}, N, length(wt.scales))
    P = zeros(Float64, N, length(wt.scales))
    for i in 1:M
        wX = transform(wt, X[:,i])
        wY = transform(wt, Y[:,i])
        C[:,:] += wX.wX .* conj.(wY.wX)
        P[:,:] += map(abs, wX.wX) .* map(abs, wX.wX)
    end
    C ./= P

    threshold = NaN;
    if nsurrogate > 0
        empC = Array{Float64, 1}[]
        # generate nsurrogate surrogate datasets and obtain coherence
        for i in 1:nsurrogate
            Cs = coherence(wt, surrogate(X); nsurrogate=0).coh
            empC = vcat(empC, map(abs, reshape(Cs, (prod(size(Cs)),)) ) )
        end
        # get threshold corresponding to the given α-level
        sort!(empC)
        i = Int32(floor((1.0-α)*length(empC)))
        threshold = empC[(i<1) ? 1 : i]
    end
    ContinuousWaveletCoherence(C, wt.scales, wt.wavelet, threshold)
end

Base.size(A::ContinuousWaveletCoherence) = size(A.coh)
Base.getindex(A::ContinuousWaveletCoherence, I::Vararg{Int, N}) where {N} = getindex(A.coh, I...);
Base.setindex!(A::ContinuousWaveletCoherence, v, I::Vararg{Int, N}) where {N} = setindex!(A.coh, v, I...)

@doc raw"""
    contourf(A::ContinuousWaveletCoherence;  drawvalid=true, shadevalid=true, kw...)

Implements a plotting helper function to plot wavelet-coherence analyses using a filled contour
plot.

If a significance test was performed using surrogate data during the analyses, the
point-wise significant areas are shown using a single (thicker) contour line.

If the optional keyword argument `drawvalid` is `true` (default), the area of valid wavelet
coefficients is marked by a thick black line. If the optional keyword arguement `shadevalid` is
also `true` (default) the invalid wavelet coefficients are shaded gray. A valid wavelet coefficient
is one, where the accociated wavelet does not overlapp with the edges of the time-series.

Additional keyword arguments are passed to the default implementation of `Plots.contourf()`.
"""
function Plots.contourf(A::ContinuousWaveletCoherence;  drawvalid=true, shadevalid=true, kw...)
    N,M = size(A)
    X = 1:N; Y = A.scales
    Z = transpose(map(abs, A.coh))
    cn = contourf(X, Y, Z;
                  linewidth=0, linetype=:heatmap,
                  xlims=(1,N), ylims=(minimum(A.scales),maximum(A.scales)),
                  kw...);
    # draw point-wise significance level if set
    if isfinite(A.α)
        contour!(X, Y, Z; levels=[A.α], linewidth=3, linecolor=:green)
    end
    # draw valid range
    if drawvalid
        x = [1.0]; y = Float64[A.scales[1]]
        for i in 1:M
            tc = cutoff_time(A.wavelet)*A.scales[i]
            if (tc < N÷2)
                push!(y, A.scales[i])
                push!(x, tc)
            end
        end
        append!(x, reverse(N.-x))
        append!(y, reverse(y))
        if shadevalid
            plot!(x,y;
                linecolor=:black, linewidth=3, label="",
                fillrange=A.scales[end], fillalpha=0.5, fillcolor=:black)
        else
            plot!(x,y;
                linecolor=:black, linewidth=3, label="")
        end
    end
    return cn
end
