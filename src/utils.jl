## Utility functions for HyperCubicRoots.jl

poor_conditioning(v, compare_tol) = any(abs.(v) .> compare_tol) ? true : false

function leading_coeff_approx_zero(polycoeff::Vector{T}; leading_tol::Float64=CUBIC_ATOL) where {T<:Real}
    if isapprox(polycoeff[end], zero(T), atol=leading_tol)
        @warn "The leading cubic coefficient is approximately zero."
        return true
    end

    return false
end

function large_scaled_coeff(polycoeff::Vector{T}; coeff_tol::Float64=1.0/CUBIC_ATOL) where {T<:Real}

    scaled_coeff = polycoeff[1:(end-1)]./polycoeff[end]
    
    if any(isinf.(scaled_coeff)) || poor_conditioning(scaled_coeff, coeff_tol)
        @warn "Large coefficients in scaled polynomial. Poor conditioning may occur."
        return true
    end

    return false
end

function test_identity(f1,f2,test_array;tol=10*eps())
    A = f1.(test_array)
    B = f2.(test_array)
    D = abs.(A .- B)
    M = maximum(D)

    return (M < tol), M, argmax(D)
end