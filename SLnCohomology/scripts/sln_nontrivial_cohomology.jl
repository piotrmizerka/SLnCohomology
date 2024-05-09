# These parameters are subject to appropriate change.
# Available options of (N,p): (3,3), (4,2).
const N = 4
const p = 2

# include("./differentials_computation/sln_laplacians.jl"); # uncomment if serialized Laplacians not available.
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using AbstractAlgebra
using BlockArrays
using Groups
using LowCohomologySOS
using JSON
using Permutations
using Serialization
using SLnCohomology
using SparseArrays

sln_laplacian_data = deserialize(joinpath(@__DIR__, "./differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
Δ = sln_laplacian_data["laplacians"]
if N == 4 # don't consider 1st cohomology for SL(4,Z) since one of the cells was not simplicial
    delete!(Δ,8)
end

@info "Laplacian loaded"

# Compute orthogonal representations without invariant vectors of a subgroup H of SL(N,p).
subgroup_rep, deg = SLnCohomology.flip_permutation_representation(N, p)
H = collect(keys(subgroup_rep))
no_inv_subspace = SLnCohomology.no_inv_subspace(H,subgroup_rep)
@assert length(no_inv_subspace) == deg # this means that we have no invariant vectors

# Induce the repesentation of H to SL(N,p)
sl_n_p_matrices = SLnCohomology.sl_n_p(N,p)
support = SLnCohomology.laplacians_support(Δ,p)
coset_data = SLnCohomology.coset_data(H,sl_n_p_matrices,p)
π = SLnCohomology.ind_rep_dict(support,subgroup_rep,coset_data,deg,p)

# Sanity check: check if π is orthogonal.
for g in support
    @assert Matrix(π[g])^(-1) == π[g]'
end

@info "Representation computed"

# Compute Laplacians evaluated on the induced representation
πΔ = Dict()
for entry in Δ
    n = entry[1]
    πΔ[n] = vcat(
        [
            hcat([SLnCohomology.representing_matrix(Δ[n][i,j],π, p) for j in 1:size(Δ[n])[2]]...)
            for i in 1:size(Δ[n])[1]
        ]...
    )
end

@info "Rep eveluation on laplacian computed"

# Show cohomologies' ranks for the above Laplacians
for entry in Δ
    n = Int8(entry[1])
    @info "rank H^"*string(div(N*(N-1),2)+N-n-1)*" = "*string(size(πΔ[n])[1]-rank(πΔ[n]))
end
