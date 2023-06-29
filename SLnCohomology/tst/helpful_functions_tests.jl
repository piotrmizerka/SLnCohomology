@testset "averaged_rep" begin
    n = 10
    G = LowCohomologySOS.cyclic_group(n)
    RG = LowCohomologySOS.group_ring(G,n)
    a, = gens(G)
    elt_list = [(a,1),(a^2,-1),(a^4,-1)]
    averaged_rep_ = averaged_rep(elt_list,RG)
    proper_average = 1//3*RG(a)-1//3*RG(a^2)-1//3*RG(a^4)
    @test averaged_rep_ == proper_average

    sl3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
    S = gens(sl3)
    S_inv =[S; inv.(S)]
    elt_list = [(S[1],1),(S[2],-1),(S[3],-1),(S[4],1)]
    RG = LowCohomologySOS.group_ring(sl3,S_inv)
    averaged_rep_ = averaged_rep(elt_list)
    proper_average = 1//4*RG(S[1])-1//4*RG(S[2])-1//4*RG(S[3])+1//4*RG(S[4])
    @test averaged_rep_ == proper_average
end

@testset "Gaussian elimination" begin
    slN = MatrixGroups.SpecialLinearGroup{4}(Int8)
    euclidean_algorithm = EuclideanAlgorithm(slN)
    e12, e13, e14, e21, e23, e24, e31, e32, e34, e41, e42, e43 = gens(slN)
    m_e12 = [1 1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    m_e23 = [1 0 0 0; 0 1 1 0; 0 0 1 0; 0 0 0 1]
    m_e31 = [1 0 0 0; 0 1 0 0; 1 0 1 0; 0 0 0 1]
    m_e42 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 1 0 1]
    @test euclidean_algorithm.gelt_from_matrix(m_e12) == e12
    @test euclidean_algorithm.gelt_from_matrix(m_e23) == e23
    @test euclidean_algorithm.gelt_from_matrix(m_e31) == e31
    @test euclidean_algorithm.gelt_from_matrix(m_e42) == e42
    @test euclidean_algorithm.gelt_from_matrix(m_e12*m_e23) == e12*e23
    @test euclidean_algorithm.gelt_from_matrix(m_e23*m_e31) == e23*e31
    @test euclidean_algorithm.gelt_from_matrix(m_e31^2) == e31^2
    @test euclidean_algorithm.gelt_from_matrix(m_e42*m_e23) == e42*e23
    @test euclidean_algorithm.gelt_from_matrix(m_e12^3) == e12^3
    @test euclidean_algorithm.gelt_from_matrix(m_e23^2*m_e42) == e23^2*e42
    @test euclidean_algorithm.gelt_from_matrix(m_e31*m_e42*m_e23) == e31*e42*e23
    @test euclidean_algorithm.gelt_from_matrix(m_e42^2*m_e12) == e42^2*e12
end
