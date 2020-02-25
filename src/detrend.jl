@doc raw"""
    detrend!(y::AbstractVector{V}, x::AbstractVector{U}=(1:length(y)))

In-place detrends the given time-series. That is, fitting a linear function into `y` given `x` and
subtracting the fit from the data in `y`.
"""
function detrend!(y::AbstractVector{V}, x::AbstractVector{U}=(1:length(y))) where {U,V}
    (length(x) != length(y)) && error("x & y must have the same length!")
    N = length(x)

    X = Matrix(zeros(Float64, N, 2))
    for i in 1:N
        X[i,:] = [1, x[i]]
    end

    res = transpose(X) * y
    res = (transpose(X) *  X) \ res
    y[:] -= (res[1] .+ (res[2] .* x))
end

@doc raw"""
    detrend!(y::AbstractArray{V,2}, x::AbstractVector{U}=(1:length(y)))

In-place detrends the given time-series. That is, fitting a linear function into every columns of
`y` given `x` and subtracting the fit from the data in the corresponding column of `y`.
"""
function detrend!(y::AbstractArray{V,2}, x::AbstractVector{U}=(1:length(y))) where {U,V}
    N,M = size(y)
    for i in 1:M
        y[:,i] = detrend!(y[:,i])
    end
    y
end

@doc raw"""
    detrend(y::AbstractVector{V}, x::AbstractVector{U}=(1:length(y)))

Detrends the given time-series. That is, fitting a linear function into `y` given `x` and
subtracting the fit from the data in `y`. The function returns the detrended time-series.
"""
function detrend(y::AbstractVector{V}, x::AbstractVector{U}=(1:length(y))) where {U,V}
    Y = copy(y)
    detrend!(Y)
    Y
end

@doc raw"""
    detrend!(y::AbstractArray{V,2}, x::AbstractVector{U}=(1:length(y)))

Detrends the given time-series. That is, fitting a linear function into every columns of
`y` given `x` and subtracting the fit from the data in the corresponding column of `y`. The fuction
returns the detrended time-series.
"""
function detrend(y::AbstractArray{V,2}, x::AbstractVector{U}=(1:length(y))) where {U,V}
    Y = copy(y)
    detrend!(Y)
    Y
end
