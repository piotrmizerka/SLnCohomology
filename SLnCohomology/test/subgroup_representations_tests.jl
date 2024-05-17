@testset "flip_permutation_representation" begin
    # Check the homomorphism condition for our reps

    # For SL(3,3):
    subgroup_rep, deg = SLnCohomology.flip_permutation_representation(3,3)
    H = keys(subgroup_rep)
    for i in 1:10
        for j in 1:10
            a, b = rand(H), rand(H) 
            @test subgroup_rep[SLnCohomology.matrix_mod_p(a*b,3)] == subgroup_rep[a]*subgroup_rep[b]
        end
    end
    # For SL(4,2):
    subgroup_rep, deg = SLnCohomology.flip_permutation_representation(4,2)
    H = keys(subgroup_rep)
    for i in 1:20
        for j in 1:20
            a, b = rand(H), rand(H)
            @test subgroup_rep[SLnCohomology.matrix_mod_p(a*b,2)] == subgroup_rep[a]*subgroup_rep[b]
        end
    end
end

@testset "subgroup_gens_expression" begin
    # Checking the orders only

    # SL(3,3)
    n, p = 3, 3
    gens_H, gens_N = SLnCohomology.symmetric_subgroups_gens(n,p)
    H_gens_expr = SLnCohomology.subgroup_gens_expression(gens_H, p)
    H = keys(H_gens_expr)
    N = keys(SLnCohomology.subgroup_gens_expression(gens_N, p))
    @test length(H) == 36
    @test length(N) == 18
    for k in N
        @test k in H
    end

    # SL(4, 2)
    n, p = 4, 2
    gens_H, gens_N = SLnCohomology.symmetric_subgroups_gens(n,p)
    H_gens_expr = SLnCohomology.subgroup_gens_expression(gens_H, p)
    H = keys(H_gens_expr)
    N = keys(SLnCohomology.subgroup_gens_expression(gens_N, p))
    @test length(H) == 576
    @test length(N) == 96
    for k in N
        @test k in H
    end
end

@testset "symmetric_subgroups_gens" begin
    # SL(3,3)
    n, p = 3, 3
    gens_H, gens_N = SLnCohomology.symmetric_subgroups_gens(n,p)

    for s in gens_H
        @test size(s) == (n,n)
    end
    for s in gens_N
        @test size(s) == (n,n)
    end

    # sanity checks for inverses of s and t
    s, t, s_inv = gens_H
    s_s_inv = SLnCohomology.matrix_mod_p(s*s_inv,p)
    s_inv_s = SLnCohomology.matrix_mod_p(s_inv*s,p)
    t2 = SLnCohomology.matrix_mod_p(t^2,p)
    @test s_s_inv == s_inv_s == t2 == I

    # sanity checks for inverses of the generators of N
    a, b, c, a_inv, b_inv, c_inv = gens_N
    a_a_inv = SLnCohomology.matrix_mod_p(a*a_inv,p)
    a_inv_a = SLnCohomology.matrix_mod_p(a_inv*a,p)
    b_b_inv = SLnCohomology.matrix_mod_p(b*b_inv,p)
    b_inv_b = SLnCohomology.matrix_mod_p(b_inv*b,p)
    c_c_inv = SLnCohomology.matrix_mod_p(c*c_inv,p)
    c_inv_c = SLnCohomology.matrix_mod_p(c_inv*c,p)
    @test a_a_inv == a_inv_a == b_b_inv == b_inv_b == c_c_inv == c_inv_c == I

    # SL(4,2)
    n, p = 4, 2
    gens_H, gens_N = SLnCohomology.symmetric_subgroups_gens(n,p)

    for s in gens_H
        @test size(s) == (n,n)
    end
    for s in gens_N
        @test size(s) == (n,n)
    end

    # sanity checks for inverses of s and t
    s,t,s_inv,t_inv = gens_H
    s_s_inv = SLnCohomology.matrix_mod_p(s*s_inv,p)
    s_inv_s = SLnCohomology.matrix_mod_p(s_inv*s,p)
    t_t_inv = SLnCohomology.matrix_mod_p(t*t_inv,p)
    t_inv_t = SLnCohomology.matrix_mod_p(t_inv*t,p)
    @test s_s_inv == s_inv_s == t_t_inv == t_inv_t == I

    # sanity checks for inverses of the generators of N
    a,b,c,d,e,f,b_inv = gens_N
    b_b_inv = SLnCohomology.matrix_mod_p(b*b_inv,p)
    b_inv_b = SLnCohomology.matrix_mod_p(b_inv*b,p)
    a2 = SLnCohomology.matrix_mod_p(a^2,p)
    c2 = SLnCohomology.matrix_mod_p(c^2,p)
    d2 = SLnCohomology.matrix_mod_p(d^2,p)
    e2 = SLnCohomology.matrix_mod_p(e^2,p)
    f2 = SLnCohomology.matrix_mod_p(f^2,p)
    @test b_b_inv == b_inv_b == I
    @test a2 == c2 == d2 == e2 == f2 == I
end
