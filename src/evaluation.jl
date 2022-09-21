## Evaluation functions for HyperCubicRoots.jl

function horner_evalpoly(degree::Integer, polycoeff::Vector{T}, x::S) where {T<:Real, S<:Real}
    y = polycoeff[degree+1]

    for i in degree:-1:1
        y = muladd(x, y, polycoeff[i])
    end

    return y
end

horner_evalpoly(polycoeff, x) = horner_evalpoly(length(polycoeff)-1, polycoeff, x)

"""
Evalutes a polynomial using Horner's method with a running error bound.

See: Higham, Nicholas J. Accuracy and stability of numerical algorithms.
Society for industrial and applied mathematics, 2002. (Chapter 5) 
"""
function horner_evalpoly_run_err_bnd(degree::Integer, polycoeff::Vector{T}, x::S) where {T<:Real, S<:Real}

    y = polycoeff[degree+1]
    mu = abs(y)/2

    for i in degree:-1:1
        y = muladd(x, y, polycoeff[i])
        mu = muladd(abs(x), mu, abs(y))
    end

    # u = eps(T)/2 # (unit round-off)
    # mu = u*(2*mu - abs(y))
    w = eps(promote_type(S,T))
    mu = w*(mu - abs(y)/2)

    return y, mu
end

horner_evalpoly_run_err_bnd(polycoeff, x::T) where T<:Number = horner_evalpoly_run_err_bnd(length(polycoeff)-1, polycoeff, x)

horner_evalpoly_run_err_bnd(polycoeff, X::T) where T<:Vector = horner_evalpoly_run_err_bnd.((polycoeff,), X)