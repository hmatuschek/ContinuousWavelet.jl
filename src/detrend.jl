
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

function detrend!(y::AbstractArray{V,2}, x::AbstractVector{U}=(1:length(y))) where {U,V}
    N,M = size(y)
    for i in 1:M
        y[:,i] = detrend!(y[:,i])
    end
    y
end

function detrend(y::AbstractVector{V}, x::AbstractVector{U}=(1:length(y))) where {U,V}
    Y = copy(y)
    detrend!(Y)
    Y
end

function detrend(y::AbstractArray{V,2}, x::AbstractVector{U}=(1:length(y))) where {U,V}
    Y = copy(y)
    detrend!(Y)
    Y
end
