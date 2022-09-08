using HyperCubicRoots
using Test

import Polynomials
import PolynomialRoots
comp_tol = 1e-3  # comparison difference tolerance to PolynomialRoots

@testset "Function evaluation" begin
    # # Testing on a known case: 
    L(z) = log(1+z)/z
    G(z) = HgF(1,1,2,-Complex(z))
    aa = 0.5:1e-2:100;
    @test HyperCubicRoots.test_identity(L, G, aa)[1] === true

    ## Cubic-related identity:
    P(z) = ∛(1 + z) + ∛(1-z);
    H(z) = 2HgF(-1/6, 1/3, 1/2, Complex(z)^2)
    bb = 0:1e-3:0.999
    test_tol=10*eps()
    @test HyperCubicRoots.test_identity(H, P, bb, tol=test_tol)[1] === true

    # We need to ensure that |z| < 1 in the argument of HgF
    # (This is done via the Kummer transformation in the cubic roots derivation.)
    @test isreal(H(1.01)) == false
    @test isapprox(P(1.01), H(1.01)) === false

    ## Evaluate against other known methods:
    cc = [-72, -36, 0, 1]

    # Hand-coded hypergeometric method:
    p1 = -12
    q1 = -36
    Del1 = q1^2 + p1^3
    s1 = Del1/p1^3
    r_hyp = [
        2q1/p1*HgF(1/3, 2/3, 1/2, s1),
    -(2q1/3p1)*HgF(2/3, 4/3, 3/2, 1/2*(1 + sqrt(s1))),
    -(2q1/3p1)*HgF(2/3, 4/3, 3/2, 1/2*(1 - sqrt(s1))),
    ]
    sort!(r_hyp)

    # Exact equality (up to ordering)
    @test r_hyp == sort(HyperCubicRoots.solve_real_cubic_roots(cc))

    # Trigonometric method
    r_trig = [
        2*sqrt(12)*cos(  π/18),
        2*sqrt(12)*cos(13π/18),
        2*sqrt(12)*cos(25π/18),
    ]
    @test all(isapprox.(sort(r_trig), r_hyp))

    # PolynomialRoots method
    r_num = PolynomialRoots.roots(cc)
    @test all(isapprox.(sort(r_num, by=real), r_hyp))
end

@testset "Numerics" begin
    ## Leading coefficient is zero:
    aa = [-417.47, 843.75, 521.67, 0.0]
    @test isempty(HyperCubicRoots.solve_real_cubic_roots(aa))

    ## Leading coefficient is zero:
    aa4 = [-417.47, 843.75, 521.67, 23.4, 0.0]
    @test isempty(HyperCubicRoots.solve_real_cubic_roots(aa4))
    
    ## Leading coefficient is nearly zero:
    bb3 = [-417.47, 843.75, 521.67, 1e-200]
    @test_throws ErrorException HyperCubicRoots.solve_real_cubic_roots(bb3)

    ## Leading coefficient is nearly zero:
    bb4 = [-417.47, 843.75, 521.67, 23.4, 1e-200]
    @test_throws ErrorException HyperCubicRoots.solve_all_quartic_roots(bb4)

    ## Three real roots example:
    cc = [-23.88, -8.04, 81.02, 13.17]
    @test length(HyperCubicRoots.solve_real_cubic_roots(cc)) == 3
    
    ## One real root example:
    dd = [81.80, -59.66, -76.47, -64.32]
    @test length(HyperCubicRoots.solve_real_cubic_roots(dd)) == 1

    ## Integer coefficient vector:
    r = HyperCubicRoots.solve_real_cubic_roots([8, 2, -5, 1])
    @test all(isapprox.(sort(r), [-1, 2, 4])) == true
end

@testset "Quartics" begin
    for pp in [
        # Constant term is zero
        [0.0, 22.4, 388.54, -442.51, 445.24], 
        # Fail cases for symmetric roots bug
        [ 936.22, 165.74, 749.85, 229.61, -553.5], 
        [-93.12, -738.32, -477.59, 237.14, 487.27],
         # Largish error compared to PolynomialRoots in Float64 precision
        BigFloat[-729.31, 346.23, -544.86, 721.04, -0.01]
    ] 

        t = HyperCubicRoots.solve_all_quartic_roots(pp)
        r = PolynomialRoots.roots(pp)

        t_r = sort(real.(t))
        t_i = sort(imag.(t))
        P_r = sort(real.(r))
        P_i = sort(imag.(r))
        @test all(isapprox.(t_r,P_r; atol=comp_tol))
        @test all(isapprox.(t_i,P_i; atol=comp_tol))
    end
end