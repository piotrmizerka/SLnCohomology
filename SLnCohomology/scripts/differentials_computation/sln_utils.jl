using Groups
using LinearAlgebra
using LowCohomologySOS
using Pkg
using Serialization
using SLnCohomology

Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

const N = 3 # hard-coded, to change if one wants other slns
sln = MatrixGroups.SpecialLinearGroup{N}(Int8)
sln_gens = gens(sln)

# Load the data from Ben
sln_bound_stab = deserialize(joinpath(@__DIR__, "./precomputed_boundaries/sl"*string(N)*"_bound_stab.sjl"))

# Extract homology degrees
differential_degrees = sort([first(x) for x in sln_bound_stab["boundaries"]])

# Compute the stabilisers as signed subgroups of SL(n,ℤ)
m_stabs = Dict()
for n in differential_degrees
    m_stabs[n] = [
        [(SLnCohomology.gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in sln_bound_stab["stabilisers"][n]
    ]
end

# Compute the supports (i.e. half_bases) for the group rings to compute the Laplacians
# - we just add to half_basis the coset elements appearing in the differentials.
d_union = Dict(k=>[one(sln)] for k in differential_degrees[2:end])
for k in differential_degrees[2:end]
    for i in eachindex(sln_bound_stab["boundaries"][k])
        for j in eachindex(sln_bound_stab["boundaries"][k][i])
            coset = sln_bound_stab["boundaries"][k][i][j]["orbit_coset_with_orientation"]
            d_union[k] = union(
                d_union[k], 
                [SLnCohomology.gelt_from_matrix(M,sln) for (M,σ) in coset]
            )
        end
    end
end
half_basis_Δ = Dict()
min_degree = differential_degrees[1]
half_basis_Δ[min_degree] = d_union[min_degree+1]
half_basis_Δ[min_degree] = unique(
    [half_basis_Δ[min_degree];inv.(half_basis_Δ[min_degree])]
)
for k in differential_degrees[3:end]
    half_basis_Δ[k-1] = union(d_union[k-1], d_union[k])
    half_basis_Δ[k-1] = unique([half_basis_Δ[k-1];inv.(half_basis_Δ[k-1])])
end

# Compute the group rings (with standard multiplciation by convolution: (1+g)(1+h)=1+g+h+gh)
# We consider only those Laplacians with reasonable size size of half_basis for their
# group rings. For now, we allow half_bases with at most 11_500 elements, which leaves
# cohomologies of degree: 1, 2, 3, 4, and 5, of which 4 and 5 remain interesting
# as we have vanishing for 1, 2, 3 from Bader-Sauer.
homology_degrees = []
for k in differential_degrees[1:end-1]
    if length(half_basis_Δ[k]) <= 11_500
        push!(homology_degrees,k)
    end
end
RG_Δ = Dict()
for k in homology_degrees
    RG_Δ[k] = LowCohomologySOS.group_ring(
        sln, half_basis_Δ[k], star_multiplication = false
    )
end

# Compute the differential as matrices over RGs
cells_number = Dict(k=>length(sln_bound_stab["stabilisers"][k]) for k in differential_degrees)
dx = Dict([(k+1,i,j) => 0//1*zero(RG_Δ[k]) for k in homology_degrees for i in 1:cells_number[k] for j in 1:cells_number[k+1]])
for k in homology_degrees
    for j in eachindex(sln_bound_stab["boundaries"][k+1])
        for summand in sln_bound_stab["boundaries"][k+1][j]
            coset_orient = [(SLnCohomology.gelt_from_matrix(M,sln),σ) for (M,σ) in summand["orbit_coset_with_orientation"]]
            summand_in_RG = summand["sign"]*SLnCohomology.averaged_rep(coset_orient,RG_Δ[k])
            dx[k+1,summand["orbit_standard_cell"],j] += summand_in_RG
        end
    end
end
d = Dict()
for k in homology_degrees
    d[k+1] = [dx[k+1,i,j] for i in 1:cells_number[k], j in 1:cells_number[k+1]]
end
