@testset "coset_data" begin
    # testing SL(2,2) ≅ S₃ ≅ D₆ = ⟨a,b|a³=b²=1,bab=a⁻¹⟩
    p = 2
    id_ = [1 0;
           0 1]
    a = [0 1;
         1 1]
    b = [0 1;
         1 0]
    a2 = SLnCohomology.matrix_mod_p(a^2,p)
    ba = SLnCohomology.matrix_mod_p(b*a,p)
    ba2 = SLnCohomology.matrix_mod_p(b*a^2,p)
    @test a2 == [1 1;
                 1 0]
    @test ba == [1 1;
                 0 1]
    @test ba2 == [1 0;
                  1 1]
    sl_2_2_matrices = [id_,a,a2,b,ba,ba2]

    C₃ = [id_,a,a2]
    coset_data = SLnCohomology.coset_data(C₃,sl_2_2_matrices,p)
    @test length(coset_data["cosets_representatives"]) == 
          length(coset_data["cosets_representatives_indices"]) == 2
    @test coset_data["elt_coset_labels"][id_] == coset_data["elt_coset_labels"][a] == coset_data["elt_coset_labels"][a2] &&
          coset_data["elt_coset_labels"][b] == coset_data["elt_coset_labels"][ba] == coset_data["elt_coset_labels"][ba2] &&
          coset_data["elt_coset_labels"][id_] != coset_data["elt_coset_labels"][b]
    @test Set(coset_data["cosets_representatives"]) == Set([id_,b]) ||
          Set(coset_data["cosets_representatives"]) == Set([id_,ba]) ||
          Set(coset_data["cosets_representatives"]) == Set([id_,ba2]) ||
          Set(coset_data["cosets_representatives"]) == Set([a,b]) ||
          Set(coset_data["cosets_representatives"]) == Set([a,ba]) ||
          Set(coset_data["cosets_representatives"]) == Set([a,ba2]) ||
          Set(coset_data["cosets_representatives"]) == Set([a2,b]) ||
          Set(coset_data["cosets_representatives"]) == Set([a2,ba]) ||
          Set(coset_data["cosets_representatives"]) == Set([a2,ba2])
    for i in eachindex(coset_data["cosets_representatives"])
        @test coset_data["cosets_representatives_indices"][coset_data["cosets_representatives"][i]] == i
    end

    C₂ = [id_,b]
    coset_data = SLnCohomology.coset_data(C₂,sl_2_2_matrices,p)
    @test length(coset_data["cosets_representatives"]) == 
          length(coset_data["cosets_representatives_indices"]) == 3
    @test coset_data["elt_coset_labels"][id_] == coset_data["elt_coset_labels"][b] &&
          coset_data["elt_coset_labels"][a] == coset_data["elt_coset_labels"][ba2] &&
          coset_data["elt_coset_labels"][a2] == coset_data["elt_coset_labels"][ba] &&
          coset_data["elt_coset_labels"][id_] != coset_data["elt_coset_labels"][a] &&
          coset_data["elt_coset_labels"][id_] != coset_data["elt_coset_labels"][a2] &&
          coset_data["elt_coset_labels"][a] != coset_data["elt_coset_labels"][a2]
    @test Set(coset_data["cosets_representatives"]) == Set([id_,a,a2]) ||
          Set(coset_data["cosets_representatives"]) == Set([id_,a,ba]) ||
          Set(coset_data["cosets_representatives"]) == Set([id_,ba,a2]) ||
          Set(coset_data["cosets_representatives"]) == Set([id_,ba,ba2]) ||
          Set(coset_data["cosets_representatives"]) == Set([b,a,a2]) ||
          Set(coset_data["cosets_representatives"]) == Set([b,a,ba]) ||
          Set(coset_data["cosets_representatives"]) == Set([b,ba,a2]) ||
          Set(coset_data["cosets_representatives"]) == Set([b,ba,ba2])
    for i in eachindex(coset_data["cosets_representatives"])
        @test coset_data["cosets_representatives_indices"][coset_data["cosets_representatives"][i]] == i
    end
end

@testset "ind_rep_dict" begin
    # testing SL(2,2) ≅ S₃ ≅ D₆ = ⟨a,b|a³=b²=1,bab=a⁻¹⟩
    p = 2
    id_ = [1 0;0 1]
    a = [0 1;1 1]
    b = [0 1;1 0]
    a2 = SLnCohomology.matrix_mod_p(a^2,p)
    ba = SLnCohomology.matrix_mod_p(b*a,p)
    ba2 = SLnCohomology.matrix_mod_p(b*a^2,p)
    sl_2_2_matrices = [id_,a,a2,b,ba,ba2]
    π = Dict(id_ => [1],b => [-1])
    invs = Dict(id_ => id_, a => a2, a2 => a, b => b, ba => ba, ba2 => ba2)
    C₂ = [id_,b]
    coset_data = SLnCohomology.coset_data(C₂,sl_2_2_matrices,p)
    π_ind = SLnCohomology.ind_rep_dict(sl_2_2_matrices,π,coset_data,1,p)
    @test π_ind[id_] == [1 0 0;
                         0 1 0;
                         0 0 1]
    for g in sl_2_2_matrices
        @assert Matrix(π_ind[g])^(-1) == π_ind[g]'
        @assert Matrix(π_ind[g])^(-1) == π_ind[invs[g]]
        for h in sl_2_2_matrices
            gh = SLnCohomology.matrix_mod_p(g*h,p)
            @assert π_ind[g]*π_ind[h] == π_ind[gh]
        end
    end
end

@testset "no_inv_subspace" begin
    H = ["a","n","y","t","h","i","nn","g"]
    x = rand((-1,1))
    y = rand((-1,1))
    z = rand((-1,1))
    π = Dict(
        "a" => [x 0 0;0 1 0;0 0 1],
        "n" => [1 0 0;0 y 0;0 0 1],
        "y" => [1 0 0;0 1 0;0 0 z],
        "t" => [1 0 0;0 1 0;0 0 1],
        "h" => [1 0 0;0 1 0;0 0 1],
        "i" => [1 0 0;0 1 0;0 0 1],
        "nn" => [1 0 0;0 1 0;0 0 1],
        "g" => [1 0 0;0 1 0;0 0 1]
    )
    proper_val = []
    if x == -1
        push!(proper_val,1)
    end
    if y == -1
        push!(proper_val,2)
    end
    if z == -1
        push!(proper_val,3)
    end
    @test SLnCohomology.no_inv_subspace(H,π) == Set(proper_val)
end

@testset "permutations_degree" begin
    @test SLnCohomology.permutations_degree([[[1 2],[3]],[[7 12],[5]]]) == 12
    @test SLnCohomology.permutations_degree([[[1,2],[30]],[[7,12],[5]]]) == 30
    @test SLnCohomology.permutations_degree([[[13 2],[3]],[[7,12],[5]]]) == 13
end

@testset "representing_matrix" begin
    p = 2
    matrices_dict, matrices_tuples_order = SLnCohomology.slnp(3,p)
    π = SLnCohomology.regular_rep(matrices_tuples_order, matrices_dict, p)
    sl3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
    e12, e13, e21, e23, e31, e32 = S = gens(sl3)
    S_inv =[S; inv.(S)]
    half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2)
    RG = LowCohomologySOS.group_ring(sl3,half_basis)
    slnp_order_ = SLnCohomology.slnp_order(3,p)
    for s in S_inv
        rep_mat = SLnCohomology.representing_matrix(RG(s),π, p)
        @test size(rep_mat) == (slnp_order_,slnp_order_) # the size of perm matrix must be equal to SLₙ(p) order
        @test length(SparseArrays.nonzeroinds(sparse(vec(rep_mat)))) == slnp_order_ # the number of nonzero indices must equal to perm degree
    end
    # representing_matrix must preserve group ring structure:
    rep_mat(ξ) = SLnCohomology.representing_matrix(ξ,π, p)
    for i in 1:10
        i, j = rand(1:12), rand(1:12)
        ξ = RG(S_inv[i]*S_inv[j]) # test on sth not too trivial
        η = RG(S_inv[j]*S_inv[i]^(-1))
        πξ = rep_mat(ξ)
        πη = rep_mat(η)
        πξ_plus_η = rep_mat(ξ+η)
        πξ_times_η = rep_mat(ξ*η)
        @test πξ_plus_η == πξ+πη
        @test πξ_times_η == πξ*πη
    end
end

@testset "standarize_permutation" begin
    @test Set(SLnCohomology.standarize_permutation([[1,2],[3]],4)) == Set([[1,2],[3],[4]])
    @test Set(SLnCohomology.standarize_permutation([[1,2],[4]],5)) == Set([[1,2],[3],[4],[5]])
end
