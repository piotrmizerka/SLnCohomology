using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using AbstractAlgebra
using BlockArrays
using Groups
using LinearAlgebraX
using LowCohomologySOS
using JSON
using Permutations
using Serialization
using SLnCohomology
using SparseArrays

# The degree n of SLₙ(ℤ) is given as a paremeter
# Available options of (n,p): (3,3), (4,2).
n = parse(Int64, ARGS[1])
p = (n == 3 ? 3 : 2)

sln_laplacian_data = deserialize(joinpath(@__DIR__, "./laplacians_computation/precomputed_laplacians/sl"*string(n)*"_laplacians.sjl"))
Δ = sln_laplacian_data["laplacians"]

@info "Laplacians loaded"

# Compute orthogonal representations without invariant vectors of a subgroup H of SL(n,p).
subgroup_rep, deg = SLnCohomology.flip_permutation_representation(n,p)
H = keys(subgroup_rep)
no_inv_subspace = SLnCohomology.no_inv_subspace(H,subgroup_rep)
@assert length(no_inv_subspace) == deg # this means that we have no invariant vectors

# Induce the repesentation of H to SL(N,p)
sl_n_p_matrices = SLnCohomology.sl_n_p(n,p)
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
    k = entry[1]
    πΔ[k] = vcat(
        [
            hcat([SLnCohomology.representing_matrix(Δ[k][i,j],π, p) for j in 1:size(Δ[k])[2]]...)
            for i in 1:size(Δ[k])[1]
        ]...
    )
end

@info "Representation evaluation on the Laplacians computed"

# Show cohomologies' ranks for the above Laplacians. 
# We use the exact rank computation from LinearAlgebraX package: https://github.com/scheinerman/LinearAlgebraX.jl
for entry in Δ
    k = Int8(entry[1])

    # Perform an additional verification of singularity for Hⁿ⁻¹.
    # For the sake of rigorousity, this is not required, however - just a double check.
    if k == div(n*(n-1),2)
        kernel_vector = nullspacex(πΔ[k])[:,1]
        for coeff in kernel_vector
            @assert typeof(coeff) <: Rational
        end
        @assert πΔ[k]*kernel_vector == [0//1 for coeff in eachindex(kernel_vector)]
    end

    # Compute the exact reduced cohomology rank.
    @info "rank H^"*string(div(n*(n-1),2)+n-k-1)*" = "*string(size(πΔ[k])[1]-rankx(πΔ[k]))
end
