@testset "Representing arrays as SLₙ(Z) elements" begin
    slN = MatrixGroups.SpecialLinearGroup{4}(Int8)
    e12, e13, e14, e21, e23, e24, e31, e32, e34, e41, e42, e43 = gens(slN)
    m_e12 = [1 1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    m_e23 = [1 0 0 0; 0 1 1 0; 0 0 1 0; 0 0 0 1]
    m_e31 = [1 0 0 0; 0 1 0 0; 1 0 1 0; 0 0 0 1]
    m_e42 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 1 0 1]
    @test SLnCohomology.gelt_from_matrix(m_e12,slN) == e12
    @test SLnCohomology.gelt_from_matrix(m_e23,slN) == e23
    @test SLnCohomology.gelt_from_matrix(m_e31,slN) == e31
    @test SLnCohomology.gelt_from_matrix(m_e42,slN) == e42
    @test SLnCohomology.gelt_from_matrix(m_e12*m_e23,slN) == e12*e23
    @test SLnCohomology.gelt_from_matrix(m_e23*m_e31,slN) == e23*e31
    @test SLnCohomology.gelt_from_matrix(m_e31^2,slN) == e31^2
    @test SLnCohomology.gelt_from_matrix(m_e42*m_e23,slN) == e42*e23
    @test SLnCohomology.gelt_from_matrix(m_e12^3,slN) == e12^3
    @test SLnCohomology.gelt_from_matrix(m_e23^2*m_e42,slN) == e23^2*e42
    @test SLnCohomology.gelt_from_matrix(m_e31*m_e42*m_e23,slN) == e31*e42*e23
    @test SLnCohomology.gelt_from_matrix(m_e42^2*m_e12,slN) == e42^2*e12
end
