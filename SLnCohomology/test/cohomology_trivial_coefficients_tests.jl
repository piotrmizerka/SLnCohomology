# The ranks of rational cohomology of SL_n
cohomology_ranks_sl_2 = Dict(0=>1, 1=>0)
cohomology_ranks_sl_3 = Dict(0=>1, 1=>0, 2=>0, 3=>0)
cohomology_ranks_sl_4 = Dict(0=>1, 1=>0, 2=>0, 3=>1, 4=>0, 5=>0, 6=>0)
cohomology_ranks_sl_5 = Dict(0=>1, 1=>0, 2=>0, 3=>0, 4=>0, 5=>0, 6=>1, 7=>0, 8=>0, 9=>0, 10=>0)

function coh_trivial_rep(N)
    # Load Laplacian
    sln_laplacian_data = deserialize(joinpath(@__DIR__, "../scripts/differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
    Δ_all = sln_laplacian_data["laplacians"]
    Δ = Dict()
    # non-simplicial arising for N=4 or 5, need to take only part of the Laplacian
    if N != 5
        Δ = Δ_all
    end
    if N == 4 # don't consider 1st cohomology for SL(4,Z) since one of the cells was not simplicial
        delete!(Δ,8)
    end

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