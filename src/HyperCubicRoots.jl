module HyperCubicRoots

# # Implementation of the real cubic root finding method described in:
# # 92.34 The cubic equation – a new look at the irreducible case
# # Author(s): I. J. ZUCKER
# # Source: The Mathematical Gazette , July 2008, Vol. 92, No. 524 (July 2008), pp. 264-268
# # Published by: The Mathematical Association
# # Stable URL: https://www.jstor.org/stable/27821778

import HypergeometricFunctions
# # Obtain the Gauss hypergeometric function ₂F₁(a, b; c; z):
F = HypergeometricFunctions._₂F₁

# Tolerance to maintain accuracy on the evaluation of the root (p(root) ≈ 0)
const CUBIC_ATOL = 10.0*sqrt(eps(Float64))

poor_conditioning(v, compare_tol) = any(abs.(v) .> compare_tol) ? true : false


"""
Return all the real roots of a cubic equation represented as a
vector of real coefficients for terms of increasing degree.
Note that any higher order coefficients (degree >= 4) are ignored.
"""
function solve_real_cubic_roots(polycoeff::Vector{T}; 
    warn_scaling::Bool=true, 
    cubic_tol::T=10*sqrt(eps(T)),
    coeff_tol::T=1/(sqrt(eps(T))) ) where {T<:Real}

    Rts = Vector{T}() # object to return with real roots inside.

    if warn_scaling && isapprox(polycoeff[4], zero(T),atol=cubic_tol)
        @warn "The leading cubic coefficient is approximately zero. Returning empty array."
        return Rts
    end

    (e, d, c) = polycoeff[1:3]./polycoeff[4]
    # We reduce f(x) = e + d*x + c*x^2 + x^3 via x == t - c/3 to
    # the standard transformed form f'(t) = 2q + 3p*t + t^3 .

    if warn_scaling && poor_conditioning((e,d,c), coeff_tol)
        @warn "Large coefficients in scaled polynomial. Poor conditioning may occur."
    end

    p = d/3 - c^2/9
    q = (c^3/27) - (c*d/6) + e/2
    Disc = q^2 + p^3  # cubic discriminant

    if Disc >= 0  # one real root: encompasses the case where p >= 0.
        u = cbrt(-q + sqrt(Disc))
        v = cbrt(-q - sqrt(Disc))

        r1 = u + v            # transformed real cubic root
        push!(Rts, r1 - c/3)  # original cubic root
    elseif Disc < 0  # three real roots : we necessarily have p < 0, hence s > 0
        h = q/p
        s = Disc/p^3
        w = sqrt(s)

        r_hyp = [
                   2*h*F(1//3, 2//3, 1//2, s),
             (-2//3)*h*F(2//3, 4//3, 3//2, 1//2*(1 + w)),
             (-2//3)*h*F(2//3, 4//3, 3//2, 1//2*(1 - w)),
        ]

        append!(Rts, r_hyp .- c/3)
    else
        error("""
        The computed discriminant is not comparable to 0;
        it may have a value of NaN or Inf due to numerical overflow.
        """)
    end

    return Rts
end

end
