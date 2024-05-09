@testset "flip_permutation_representation" begin
    # Linear rep of order 2 of Q₈ embedded in SL(3,3).
    # Elements of Q₈ grouped by conjugacy classes:
    Q8 = [
        [1 0 0;
         0 1 0;
         0 0 1], # (order 1),
        [2 0 0;
         0 2 0;
         0 0 1], # (order 2),
        [2 2 0;
         2 1 0;
         0 0 1], 
        [1 1 0;
         1 2 0;
         0 0 1], # (order 4),
        [0 2 0;
         1 0 0;
         0 0 1], 
        [0 1 0;
         2 0 0;
         0 0 1], # (order 4),
        [1 2 0;
         2 2 0;
         0 0 1], 
        [2 1 0;
         1 1 0;
         0 0 1] # (order 4)
    ]
    function order_two_q8(g)
        if g == [0 2 0;1 0 0;0 0 1] || g == [0 1 0;2 0 0;0 0 1] || 
           g == [2 2 0;2 1 0;0 0 1] || g == [1 1 0;1 2 0;0 0 1]
            return reshape([-1],1,1)
        else
            return reshape([1],1,1)
        end
    end
    subgroup_rep_proper = Dict()
    for g in Q8
        subgroup_rep_proper[g] = order_two_q8(g)
    end

    subgroup_rep, deg = SLnCohomology.flip_permutation_representation(3,3)
    
    @test (subgroup_rep, deg) == (subgroup_rep_proper, 1)

    # Check the homomorphism condition for our reps
    # For SL(3,3):
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
    N, p = 3, 3
    gens_H, gens_index_two_H = SLnCohomology.symmetric_subgroups_gens(N,p)
    H_gens_expr = SLnCohomology.subgroup_gens_expression(gens_H, p)
    H = keys(H_gens_expr)
    index_two_H = keys(SLnCohomology.subgroup_gens_expression(gens_index_two_H, p))
    @test length(H) == 8
    @test length(index_two_H) == 4
    for k in index_two_H
        @test k in H
    end

    # SL(4, 2)
    N, p = 4, 2
    gens_H, gens_index_two_H = SLnCohomology.symmetric_subgroups_gens(N,p)
    H_gens_expr = SLnCohomology.subgroup_gens_expression(gens_H, p)
    H = keys(H_gens_expr)
    index_two_H = keys(SLnCohomology.subgroup_gens_expression(gens_index_two_H, p))
    @test length(H) == 576
    @test length(index_two_H) == 288
    for k in index_two_H
        @test k in H
    end
end

@testset "symmetric_subgroups_gens" begin
    # SL(3,3)
    N, p = 3, 3
    gens_H, gens_index_two_H = SLnCohomology.symmetric_subgroups_gens(N,p)

    for s in gens_H
        @test size(s) == (N,N)
    end
    for s in gens_index_two_H
        @test size(s) == (N,N)
    end

    # sanity checks for inverses of s and t
    a, b, a_inv, b_inv = gens_H
    a_a_inv = SLnCohomology.matrix_mod_p(a*a_inv,p)
    a_inv_a = SLnCohomology.matrix_mod_p(a_inv*a,p)
    b_b_inv = SLnCohomology.matrix_mod_p(b*b_inv,p)
    b_inv_b = SLnCohomology.matrix_mod_p(b_inv*b,p)
    @test a_a_inv == a_inv_a == b_b_inv == b_inv_b == I

    # sanity checks for inverses of the generators of index_two_H
    c, c_inv = gens_index_two_H
    c_c_inv = SLnCohomology.matrix_mod_p(c*c_inv,p)
    c_inv_c = SLnCohomology.matrix_mod_p(c_inv*c,p)
    @test c_c_inv == c_inv_c == I

    # SL(4,2)
    N, p = 4, 2
    gens_H, gens_index_two_H = SLnCohomology.symmetric_subgroups_gens(N,p)

    for s in gens_H
        @test size(s) == (N,N)
    end
    for s in gens_index_two_H
        @test size(s) == (N,N)
    end

    # sanity checks for inverses of s and t
    s,t,s_inv,t_inv = gens_H
    s_s_inv = SLnCohomology.matrix_mod_p(s*s_inv,p)
    s_inv_s = SLnCohomology.matrix_mod_p(s_inv*s,p)
    t_t_inv = SLnCohomology.matrix_mod_p(t*t_inv,p)
    t_inv_t = SLnCohomology.matrix_mod_p(t_inv*t,p)
    @test s_s_inv == s_inv_s == t_t_inv == t_inv_t == I

    # sanity checks for inverses of the generators of index_two_H
    a,b,c,d,e,f,g,b_inv,c_inv,g_inv = gens_index_two_H
    b_b_inv = SLnCohomology.matrix_mod_p(b*b_inv,p)
    b_inv_b = SLnCohomology.matrix_mod_p(b_inv*b,p)
    c_c_inv = SLnCohomology.matrix_mod_p(c*c_inv,p)
    c_inv_c = SLnCohomology.matrix_mod_p(c_inv*c,p)
    g_g_inv = SLnCohomology.matrix_mod_p(g*g_inv,p)
    g_inv_g = SLnCohomology.matrix_mod_p(g_inv*g,p)
    a2 = SLnCohomology.matrix_mod_p(a^2,p)
    d2 = SLnCohomology.matrix_mod_p(d^2,p)
    e2 = SLnCohomology.matrix_mod_p(e^2,p)
    f2 = SLnCohomology.matrix_mod_p(f^2,p)
    @test b_b_inv == b_inv_b == c_c_inv == c_inv_c == g_g_inv == g_inv_g == I
    @test a2 == d2 == e2 == f2 == I
end
