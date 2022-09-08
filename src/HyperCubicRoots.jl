module HyperCubicRoots

# # Obtain the Gauss hypergeometric function ₂F₁(a, b; c; z):
import HypergeometricFunctions
HgF = HypergeometricFunctions._₂F₁

export HgF,
       solve_real_cubic_roots, solve_all_cubic_roots,
       solve_quadratic_roots,
       solve_all_quartic_roots

# Tolerance to maintain accuracy on the evaluation of the root (p(root) ≈ 0)
const CUBIC_ATOL = 10.0*sqrt(eps(Float64))
const QUARTIC_ATOL = 10.0*sqrt(eps(Float64))

include("utils.jl")
include("quadratic.jl")
include("cubic.jl")
include("quartic.jl")

include("evaluation.jl")

end
