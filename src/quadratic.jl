# # quadratic.jl
# # Source: PolynomialRoots.jl: solve_quadratic_eq

function solve_quadratic_roots(polycoeff::Vector{T}) where {T<:Real}

    Rts = Vector{Complex{T}}()

    if length(polycoeff) < 3
        @error "There are too few coefficients given."
        return Rts
    end

    (c, b, a) = polycoeff

    if isapprox(a,zero(T))
        @warn "The leading cubic coefficient is approximately zero. Returning empty array."
        return Rts
    end

    # This is copied from PolynomialRoots.jl:
    Δ = sqrt(Complex(b*b - 4*a*c))
    if real(conj(b)*Δ) >= 0
        x0 = -(b + Δ)/2
    else
        x0 = -(b - Δ)/2
    end
    if x0 == 0
        x1 = x0
    else
        x1 = c / x0 # Viete's formula
        x0 = x0 / a
    end
    # end copied section.
    push!(Rts,x0,x1)
    return Rts
end