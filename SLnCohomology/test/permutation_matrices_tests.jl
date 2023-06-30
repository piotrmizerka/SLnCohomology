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

@testset "slnp" begin
    N = 4
    p = 2
    matrices_dict, matrices_tuples_order = SLnCohomology.slnp(N,p)
    @test length(matrices_dict) == length(matrices_tuples_order) == SLnCohomology.slnp_order(N,p)
    for i in eachindex(matrices_tuples_order)
        @test matrices_dict[matrices_tuples_order[i]] == i
    end
end

@testset "projection" begin
    @test SLnCohomology.projection([12 7;-10 5],3) == [0 1;2 2]
    @test SLnCohomology.projection([-12 -7;10 -5],3) == [0 2;1 1]
    @test SLnCohomology.projection([101 70;10 51],5) == [1 0;0 1]
    @test SLnCohomology.projection([-11 45;23 67],7) == [3 3;2 4]
    @test SLnCohomology.projection([8 7 4;1 50 8;-1 -2 -30],5) == [3 2 4;1 0 3;4 3 0]
end

@testset "permutation_matrix" begin
    identity_perm = PermutationGroups.Perm([1 2 3 4])
    perm1 = PermutationGroups.Perm([3 2 1])
    perm2 = PermutationGroups.Perm([1 3 2 4])
    perm3 = PermutationGroups.Perm([4 3 1 2])
    perm4 = PermutationGroups.Perm([4 1 5 3 2])
    @test SLnCohomology.permutation_matrix(identity_perm) == [
        0 0 1;
        0 1 0;
        1 0 0
    ]
    @test SLnCohomology.permutation_matrix(perm1) == [
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1
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

@testset "representing_matrix" begin
    p = 2
    matrices_dict, matrices_tuples_order = SLnCohomology.slnp(3,p)
    perm_mats = SLnCohomology.permutation_matrices(matrices_tuples_order, matrices_dict, p)
    sl4 = MatrixGroups.SpecialLinearGroup{3}(Int8)
    e12, e13, e21, e23, e31, e32 = S = gens(sl4)
    S_inv =[S; inv.(S)]
    RG = LowCohomologySOS.group_ring(sl4,S_inv)
    slnp_order_ = SLnCohomology.slnp_order(3,p)
    for s in S_inv
        rep_mat = SLnCohomology.representing_matrix(RG(s),p,perm_mats)
        @test size(rep_mat) == (slnp_order_,slnp_order_) # the sieze of perm matrix must be equal to SLₙ(p) order
        @test length(SparseArrays.nonzeroinds(rep_mat)) == slnp_order_ # the number of nonzero indices must equal to perm degree
    end
    # representing_matrix must preserve group ring structure:
    for i in 1:10
        i, j = rand(1:24), rand(1:24)
        ξ = RG(S_inv[i]^rand(2:3)*S_inv[j]*S_inv^(-1)) # test on sth not too trivial
        η = RG(S_inv[j]^rand(2:3)*S_inv[i]*S_inv^(-2))
        πξ = SLnCohomology.representing_matrix(ξ)
        πη = SLnCohomology.representing_matrix(η)
        πξ_plus_η = SLnCohomology.representing_matrix(ξ+η)
        πξ_times_η = SLnCohomology.representing_matrix(ξ*η)
        @test πξ_plus_η == πξ+πη
        @test πξ_times_η == πξ*πη
    end
end