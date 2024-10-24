# The ranks of rational cohomology of SL_n
cohomology_ranks_sl_2 = Dict(0=>1, 1=>0)
cohomology_ranks_sl_3 = Dict(0=>1, 1=>0, 2=>0, 3=>0)
cohomology_ranks_sl_4 = Dict(0=>1, 1=>0, 2=>0, 3=>1, 4=>0, 5=>0, 6=>0)

function coh_trivial_rep(N)
    # Load Laplacian
    sln_laplacian_data = deserialize(joinpath(@__DIR__, "../scripts/laplacians_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
    Δ = sln_laplacian_data["laplacians"]

    # Compute Laplacians evaluated on the induced representation (i.e. compute cohomology)
    # Here: trivial representation, so all group elements sent to the identity matrix
    πΔ = Dict()
    for entry in Δ
        n = entry[1]
        πΔ[n] = vcat(
            [
                hcat([sum(Δ[n][i,j].coeffs) for j in 1:size(Δ[n])[2]]...)
                for i in 1:size(Δ[n])[1]
            ]...
        )
    end

    coh_ranks = Dict()
    for entry in Δ
        n = Int8(entry[1])
        degree = div(N*(N-1),2)+N-n-1
        rk = size(πΔ[n])[1]-rank(πΔ[n])
        coh_ranks[degree] = rk
    end
    return coh_ranks
end

@testset "SL_2_cohomology_trivial_coeffs" begin
    coh_ranks = coh_trivial_rep(2)
    for degree in keys(coh_ranks)
        @test coh_ranks[degree] == cohomology_ranks_sl_2[degree]
    end
end

@testset "SL_3_cohomology_trivial_coeffs" begin
    coh_ranks = coh_trivial_rep(3)
    for degree in keys(coh_ranks)
        @test coh_ranks[degree] == cohomology_ranks_sl_3[degree]
    end
end

@testset "SL_4_cohomology_trivial_coeffs" begin
    coh_ranks = coh_trivial_rep(4)
    for degree in keys(coh_ranks)
        @test coh_ranks[degree] == cohomology_ranks_sl_4[degree]
    end
end