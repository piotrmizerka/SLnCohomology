using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using LowCohomologySOS
using Serialization
using SLnCohomology

# This parameter (hard-coded) is subject to appropriate change.
N = 3

# Compute the differentials
include("sln_boundary_data.jl");

# Sanity checks for vanishing compositions of differentials
# we can do sanity check for consecutive degrees only
consecutive_differential_degrees = []
for i in 1:(length(homology_degrees)-1)
    if homology_degrees[i+1] == homology_degrees[i]+1
        push!(consecutive_differential_degrees,(homology_degrees[i]+1,homology_degrees[i]+2))
    end
end
for pair in consecutive_differential_degrees
    n = pair[2]
    d_n_minus_1 = LowCohomologySOS.embed.(identity, d[n-1], Ref(RG_Δ[n-1]))
    d_n = LowCohomologySOS.embed.(identity, d[n], Ref(RG_Δ[n-1]))
    @assert d_n_minus_1*d_n == [zero(RG_Δ[n-1]) for i in 1:cells_number[n-2],j in 1:cells_number[n]]
end

# The stabiliser parts which we have to add to get free modules.
# Hopefully the stabilisers' elements belong to the half_bases.
stab_part_dim = Dict()
for n in homology_degrees
    k_n = cells_number[n]
    stab_part_dim[n] = [
        i == j ? one(RG_Δ[n])-SLnCohomology.averaged_rep(m_stabs[n][i],RG_Δ[n]) : zero(zero(RG_Δ[n])) 
        for i in 1:k_n,j in 1:k_n
    ]
end

# Compute the Laplacians (time consuming, can take about ??).
# At the end, we embed the Laplacian into RG_star, the group ring
# with the same basis as RG but with twisted multiplciation, i.e.
# (1+g)(1+h)=1+g+h+g^(-1)h. This is needed to translate hermitian squares
# to the standard ones for the definition of the semi-definite optimization problem
# (solvers prefer standard squares to hermitian :).
# This has no effect on the shape of the Laplacian since we just embed it.
RG_Δ_star = Dict()
for n in homology_degrees
    RG_Δ_star[n] = LowCohomologySOS.group_ring(sln, half_basis_Δ[n], star_multiplication = true)
end
Δ = Dict()
for pair in consecutive_differential_degrees
    n = pair[1]
    d_n = LowCohomologySOS.embed.(identity, d[n], Ref(RG_Δ[n]))
    d_n_plus_1 = LowCohomologySOS.embed.(identity, d[n+1], Ref(RG_Δ[n]))
    Δ[n] = d_n'*d_n+d_n_plus_1*d_n_plus_1'+stab_part_dim[n]
    Δ[n] = LowCohomologySOS.embed.(identity, Δ[n], Ref(RG_Δ_star[n]))
end
if homology_degrees[1] == differential_degrees[1]
    nx = homology_degrees[1]
    d_nx_plus_1 = LowCohomologySOS.embed.(identity, d[nx+1], Ref(RG_Δ[nx]))
    Δ[nx] = d_nx_plus_1*d_nx_plus_1'+stab_part_dim[nx]
end
if homology_degrees[end] == differential_degrees[end]-1
    nx = differential_degrees[end]
    d_nx = LowCohomologySOS.embed.(identity, d[nx], Ref(RG_Δ[nx-1]))
    k_nx = cells_number[nx]
    stab_part_dim[nx] = [
        i == j ? one(RG_Δ[nx-1])-SLnCohomology.averaged_rep(
            m_stabs[nx][i],RG_Δ[nx-1]) : zero(zero(RG_Δ[nx-1])
            ) 
        for i in 1:k_nx,j in 1:k_nx
    ]
    Δ[nx] = d_nx'*d_nx+stab_part_dim[nx]
    push!(homology_degrees,nx)
end

# A sanity check: check if the Laplacians are hermitian.
for entry in Δ
    @assert entry[2]' == entry[2]
end

laplacian_data = Dict()
laplacian_data["laplacians"] = Δ
laplacian_data["differential_degrees"] = differential_degrees
serialize(joinpath(@__DIR__, "./precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"), laplacian_data)
