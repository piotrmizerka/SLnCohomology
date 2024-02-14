@testset "averaged_rep" begin
    n = 10
    G = cyclic_group(n)
    RG = LowCohomologySOS.group_ring(G,n)
    a, = gens(G)
    elt_list = [(a,1),(a^2,-1),(a^4,-1)]
    averaged_rep_ = SLnCohomology.averaged_rep(elt_list,RG)
    proper_average = 1//3*RG(a)-1//3*RG(a^2)-1//3*RG(a^4)
    @test averaged_rep_ == proper_average

    sl3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
    S = gens(sl3)
    S_inv =[S; inv.(S)]
    elt_list = [(S[1],1),(S[2],-1),(S[3],-1),(S[4],1)]
    RG = LowCohomologySOS.group_ring(sl3,S_inv)
    averaged_rep_ = SLnCohomology.averaged_rep(elt_list,RG)
    proper_average = 1//4*RG(S[1])-1//4*RG(S[2])-1//4*RG(S[3])+1//4*RG(S[4])
    @test averaged_rep_ == proper_average
end

@testset "matrix_mod_p" begin
    @test SLnCohomology.matrix_mod_p([12 7;-10 5],3) == [0 1;2 2]
    @test SLnCohomology.matrix_mod_p([-12 -7;10 -5],3) == [0 2;1 1]
    @test SLnCohomology.matrix_mod_p([101 70;10 51],5) == [1 0;0 1]
    @test SLnCohomology.matrix_mod_p([-11 45;23 67],7) == [3 3;2 4]
    @test SLnCohomology.matrix_mod_p([8 7 4;1 50 8;-1 -2 -30],5) == [3 2 4;1 0 3;4 3 0]
end

@testset "permutation_matrix" begin
    identity_perm = PermutationGroups.Perm([1,2,3,4])
    perm1 = PermutationGroups.Perm([3,2,1])
    perm2 = PermutationGroups.Perm([1,3,2,4])
    perm3 = PermutationGroups.Perm([4,3,1,2])
    perm4 = PermutationGroups.Perm([4,1,5,3,2])
    @test SLnCohomology.permutation_matrix(identity_perm) == [
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1    
    ]
    @test SLnCohomology.permutation_matrix(perm1) == [
        0 0 1;
        0 1 0;
        1 0 0
    ]
    @test SLnCohomology.permutation_matrix(perm2) == [
        1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1
    ]
    @test SLnCohomology.permutation_matrix(perm3) == [
        0 0 1 0;
        0 0 0 1;
        0 1 0 0;
        1 0 0 0
    ]
    @test SLnCohomology.permutation_matrix(perm4) == [
        0 1 0 0 0;
        0 0 0 0 1;
        0 0 0 1 0;
        1 0 0 0 0;
        0 0 1 0 0
    ]
end

# This test takes up to 15 min on a standard laptop
@testset "regular_rep" begin
    N = 3
    p = 3

    slnp_dict, matrices_fixed_order = SLnCohomology.slnp(N,p)
    π_reg = SLnCohomology.regular_rep(matrices_fixed_order, slnp_dict, p)

    sln_laplacian_data = deserialize(joinpath(@__DIR__, "../scripts/differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
    Δ = sln_laplacian_data["laplacians"]

    π_reg_Δ = Dict()
    for entry in Δ
        n = entry[1]
        π_reg_Δ[n] = vcat(
            [
                hcat([SLnCohomology.representing_matrix(Δ[n][i,j],π_reg, p) for j in 1:size(Δ[n])[2]]...)
                for i in 1:size(Δ[n])[1]
            ]...
        )
    end

    # A sanity check - since H^1(SL(3,Z),π) = 0, we shall get the full rank for the perm repr of Δ₄:
    @test π_reg_Δ[4]' == π_reg_Δ[4]
    @test size(π_reg_Δ[4])[1]-rank(π_reg_Δ[4]) == 0 # rank(H^1) == 0

    # A sanity check: since the permutation representation has a nontrivial fixed pont set, 
    # we shall get nontrivial zero cohomology for this representation.
    @test π_reg_Δ[5]' == π_reg_Δ[5]
    @test size(π_reg_Δ[5])[1]-rank(π_reg_Δ[5]) > 0 # non-full rank - means nontrivial 0-cohomology

    # H^2 - we get nontrivial cohomology!
    @assert π_reg_Δ[3]' == π_reg_Δ[3]
    @test size(π_reg_Δ[3])[1]-rank(π_reg_Δ[3]) > 0 # non-full rank - means nontrivial 2-cohomology
end

@testset "slnp" begin
    N = 4
    p = 2
    matrices_dict, matrices_tuples_order = SLnCohomology.slnp(N,p)
    @test length(matrices_dict) == length(matrices_tuples_order) == SLnCohomology.slnp_order(N,p)
    for i in eachindex(matrices_tuples_order)
        @test matrices_dict[matrices_tuples_order[i]] == i
    end
end

@testset "slnp_order" begin
    @test SLnCohomology.slnp_order(2,2) == 6
    @test SLnCohomology.slnp_order(2,3) == 24
    @test SLnCohomology.slnp_order(2,5) == 120
    @test SLnCohomology.slnp_order(2,7) == 336
    @test SLnCohomology.slnp_order(2,11) == 1320
    @test SLnCohomology.slnp_order(3,2) == 168
    @test SLnCohomology.slnp_order(3,3) == 5616
    @test SLnCohomology.slnp_order(3,5) == 372_000
    @test SLnCohomology.slnp_order(3,7) == 5_630_688
    @test SLnCohomology.slnp_order(4,2) == 20_160
    @test SLnCohomology.slnp_order(4,3) == 12_130_560
    @test SLnCohomology.slnp_order(4,5) == 29_016_000_000
    @test SLnCohomology.slnp_order(5,2) == 9_999_360
    @test SLnCohomology.slnp_order(5,3) == 237_783_237_120
end
