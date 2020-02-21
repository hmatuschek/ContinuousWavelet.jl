using FFTW: r2r, r2r!, R2HC, HC2R
using Random


function surrogate(A::AbstractArray{Float64, 1})
    N = length(A)
    # compute real -> half-complex DFT
    res = r2r(A, R2HC, 1)
    # shuffle phases in half-complex representation
    ϕ = rand(Float64, N÷2)
    for i in 1:(N÷2)
        z = res[1+i] + 1im*res[end-i];
        z = abs(z)*(exp(2im*π*ϕ[i]))
        res[1+i] = real(z)
        res[end-i] = imag(z)
    end
    # compute inverse half-complex -> real DFT inplace
    r2r!(res, HC2R, 1)
    return res;
end

function surrogate!(A::AbstractArray{Float64, 1})
    N = length(A)
    # compute real -> half-complex DFT
    r2r!(A, R2HC, 1)
    # shuffle phases in half-complex representation
    ϕ = rand(Float64, N÷2)
    for i in 1:(N÷2)
        z = A[1+i] + 1im*A[end-i];
        z = abs(z)*(exp(2im*π*ϕ[i]))
        A[1+i] = real(z)
        A[end-i] = imag(z)
    end
    # compute inverse half-complex -> real DFT inplace
    r2r!(A, HC2R, 1)
end


function surrogate(A::AbstractArray{Float64, 2})
    N,M = size(A)
    # compute real -> half-complex DFT
    res = r2r(A, R2HC, 1)
    # shuffle phases in half-complex representation
    ϕ = rand(Float64, M)
    for i in 1:(N÷2)
        z = res[1+i,:] + 1im*res[end-i,:];
        z = map(abs,z).*(exp.(2im*π*ϕ))
        res[1+i,:] = map(real,z)
        res[end-i,:] = map(imag,z)
    end
    # compute inverse half-complex -> real DFT inplace
    r2r!(res, HC2R, 1)
    return res;
end

function surrogate!(A::AbstractArray{Float64, 2})
    N,M = size(A)
    # compute real -> half-complex DFT inplace
    r2r!(A, R2HC, 1)
    # shuffle phases in half-complex representation
    ϕ = rand(Float64, M)
    for i in 1:(N÷2)
        z = A[1+i,:] + 1im*A[end-i,:];
        z = map(abs,z)*(exp.(2im*π*ϕ))
        A[1+i,:] = map(real,z)
        A[end-i,:] = map(imag,z)
    end
    # compute inverse half-complex -> real DFT inplace
    r2r!(A, HC2R, 1)
    return A;
end
