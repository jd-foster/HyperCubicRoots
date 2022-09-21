# # Implementation of the quartic root finding method described in:
# # A Geometric Interpretation of the Solution of the General Quartic Polynomial
# # Author(s): William M. Faucette
# # Source: The American Mathematical Monthly , Jan., 1996, Vol. 103, No. 1 (Jan., 1996), pp. 51-57
# # Published by: The Mathematical Association of America
# # Stable URL: https://www.jstor.org/stable/2975214

const S = [-1 1 1; -1 -1 -1; 1 1 -1; 1 -1 1] # permuting signs on the roots

"""
Defined by the taking the positive square root relationship
    (-x)^(1/2) = xsign(x)*sqrt(abs(x))
"""
function xsign(x::R) where {R<:Real}

    if sign(x) == 1.0
        return im
    elseif sign(x) == -1.0
        return 1.0
    else
        error("Sign of input value not comparable to ±1.0")
        return nothing
    end

end

"""
Return all the roots of a quartic equation, represented as a
vector of real coefficients for terms of increasing degree.
Note that any higher order coefficients (degree >= 5) are ignored.

Theory described in Faucette (1996).
"""
function solve_all_quartic_roots(polycoeff::Vector{T}; realroots::Bool=false) where {T<:Real}

    Rts = realroots ? Vector{float(T)}() : Vector{Complex{float(T)}}() # object to return with roots inside.

    if length(polycoeff) < 5
        @error "Too few polynomial coefficients were given as a quartic."
        return Rts
    end

    if leading_coeff_approx_zero(polycoeff, leading_tol=eps(Float64(0.0)))
        return Rts
    end

    a = polycoeff[1:4]./polycoeff[5]

    # We reduce p(z) = a1 + a2*z + a3*z^2 + a4*z^3 + z^4 via z == x - a4/4 to
    # the standard transformed form P(x) = r + q*x + p*x^2 + x^4

    p = -3*(a[4]^2)/8   + a[3]
    q =    (a[4]^3)/8   - a[4]*a[3]/2 + a[2]
    r = -3*(a[4]^4)/256 + (a[4]^2)*a[3]/16 - a[4]*a[2]/4 + a[1]

    # The "resolvent cubic of P(x)" is h(z) = z^3 - 2p*z^2 + (p^2 - 4r)z + q^2.
    # We find it's roots:
    resolvent_coeff = [q^2, p^2 - 4r, -2p, 1]
    resolvent_roots = solve_all_cubic_roots(resolvent_coeff)

    alpha = convert(float(T),resolvent_roots[1])  # ensure first root is Real.
    beta, gamma = resolvent_roots[2:3]
    
    if !isreal(beta) && !isreal(gamma) && sign(alpha) == -1.0
        A = sqrt(-alpha)
        B = sqrt(-beta)
        C = sign(q)*sqrt(-gamma)
    elseif isreal(beta) && isreal(gamma)
        A = xsign(alpha)*sqrt(abs(alpha))
        B = xsign(real(beta))*sqrt(abs(real(beta)))
        C_sign = sign(q)/(xsign(alpha)*xsign(real(beta)))
        C = C_sign*sqrt(abs(real(gamma)))
    else
        error("Unrecognised case.")
    end

    R_res = S*([A; B; C])/2
    # R_res = S*(sqrt.(-resolvent_roots))/2

    # # A, B, C = sqrt.(-resolvent_roots)
    # # x1 = (-A + B + C)/2
    # # x2 = (-A - B - C)/2
    # # x3 = ( A + B - C)/2
    # # x4 = ( A - B + C)/2
    
    Rts = R_res .- a[4]/4

    if realroots && !all(isreal(Rts))
        @warn "The `realroots`` flag is `true`, but complex roots were found. Returning complex roots."
        return convert(Vector{Complex{float(T)}}, Rts)
    end

    return Rts
end
