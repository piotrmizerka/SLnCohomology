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
        @test Matrix(π_ind[g])^(-1) == π_ind[g]' == π_ind[invs[g]]
        for h in sl_2_2_matrices
            gh = SLnCohomology.matrix_mod_p(g*h,p)
            @test π_ind[g]*π_ind[h] == π_ind[gh]
        end
    end
end
