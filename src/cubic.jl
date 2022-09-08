# # Implementation of the real cubic root finding method described in:
# # 92.34 The cubic equation – a new look at the irreducible case
# # Author(s): I. J. ZUCKER
# # Source: The Mathematical Gazette , July 2008, Vol. 92, No. 524 (July 2008), pp. 264-268
# # Published by: The Mathematical Association of America
# # Stable URL: https://www.jstor.org/stable/27821778


"""
Return all the real roots of a cubic equation represented as a
vector of real coefficients for terms of increasing degree.
Note that any higher order coefficients (degree >= 4) are ignored.

Implementation of the real cubic root finding method of Zucker (2008).
"""
function solve_real_cubic_roots(polycoeff::Vector{T}) where {T<:Real}

    Rts = Vector{float(T)}() # object to return with real roots inside.

    if length(polycoeff) < 4
        @error "There are too few coefficients given."
        return Rts
    end

    if leading_coeff_approx_zero(polycoeff; leading_tol=eps(Float64(0.0)))
        return Rts
    end

    (e, d, c) = polycoeff[1:3]./polycoeff[4]
   
    # We reduce f(x) = e + d*x + c*x^2 + x^3 via x == t - c/3 to
    # the standard transformed form f'(t) = 2q + 3p*t + t^3 :
    p = d/3 - c^2/9
    q = (c^3/27) - (c*d/6) + e/2
    Disc = q^2 + p^3  # cubic discriminant

    if Disc >= 0  # one real root: encompasses the case where p >= 0.
        u = cbrt(-q + sqrt(Disc))
        v = cbrt(-q - sqrt(Disc))

        r1 = u + v            # transformed real cubic root
        push!(Rts, r1 - c/3)  # original cubic root
    elseif Disc < 0  # three real roots : we necessarily have p < 0, hence 0 < s < 1
        h = q/p
        s = Disc/p^3
        w = sqrt(s)

        r_hyp = [
                   2*h*HgF(1//3, 2//3, 1//2, s),
             (-2//3)*h*HgF(2//3, 4//3, 3//2, 1//2*(1 + w)),
             (-2//3)*h*HgF(2//3, 4//3, 3//2, 1//2*(1 - w)),
        ]

        append!(Rts, r_hyp .- c/3)
    else
        error("""
        The computed discriminant $(Disc) is not comparable to 0;
        it may have a value of NaN or Inf due to numerical issues (e.g. overflow).
        """)
    end

    return Rts
end

function solve_all_cubic_roots(polycoeff::Vector{T}) where {T<:Real}

    Rts = Vector{Complex{Float64}}()

    real_roots = solve_real_cubic_roots(polycoeff)
    append!(Rts, real_roots)

    if length(real_roots) == 3
        return Rts
    elseif length(real_roots) == 1
        r = real_roots[1]
        c0, c1, c2 = polycoeff[1:3]./polycoeff[4]
        # If we know the one real root r of c0 + c1*z + c2*z^2 + z^3 = 0,
        # then we can solve a quadratic -c0/r + (c2 + r)*x + x^2 = 0 (or multiply through by r)
        # to obtain s1 and s2, the complex conjugate roots, since
        # c0 = -r*s1*s2 and c2 = -(r + s1 + s2), 
        # and we obtain the sum and product of s1 and s2:
        complex_roots = solve_quadratic_roots([-c0, r*(r + c2), r]) # return type Vector{Complex{Float64}}
        Rts = append!(Rts, complex_roots)
        return Rts
    end

    ## Falls through to:
    error("Incorrect number of cubic roots returned: $(Rts).")
    return Rts
    
end