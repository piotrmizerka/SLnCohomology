using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using LowCohomologySOS
using Serialization
using SLnCohomology

# The degree of SLₙ(ℤ)
n = parse(Int64, ARGS[1])

# The boundary and stabiliser data
cells_sln = SLnCohomology.cells_sln(n)
oriented_cells_sln = SLnCohomology.oriented_cells_dict(cells_sln)
boundaries = SLnCohomology.boundaries_dict(oriented_cells_sln)
stabilisers = SLnCohomology.stabilisers_dict(oriented_cells_sln)

# Extract homology degrees
differential_degrees = sort([first(x) for x in boundaries])

# Compute the stabilisers as signed subgroups of SL(n,ℤ)
m_stabs = Dict()
sln = MatrixGroups.SpecialLinearGroup{n}(Int8)
for k in differential_degrees
    m_stabs[k] = [
        [(SLnCohomology.gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in stabilisers[k]
    ]
end

# Compute the supports (i.e. half_bases) for the group rings to compute the Laplacians
# - we just add to half_basis the coset elements appearing in the differentials.
d_union = Dict(k=>[one(sln)] for k in differential_degrees[2:end])
for k in differential_degrees[2:end]
    for i in eachindex(boundaries[k])
        for j in eachindex(boundaries[k][i])
            coset = boundaries[k][i][j]["orbit_coset_with_orientation"]
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

# Compute the differential as matrices over RGs.
# We store them as matrices over Rationals to perform exact operations.
cells_number = Dict(k=>length(stabilisers[k]) for k in differential_degrees)
dx = Dict([(k+1,i,j) => 0//1*zero(RG_Δ[k]) for k in homology_degrees for i in 1:cells_number[k] for j in 1:cells_number[k+1]])
for k in homology_degrees
    for j in eachindex(boundaries[k+1])
        for summand in boundaries[k+1][j]
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

# Sanity checks for vanishing compositions of differentials.
consecutive_differential_degrees = []
for i in 1:(length(homology_degrees)-1)
    if homology_degrees[i+1] == homology_degrees[i]+1
        push!(consecutive_differential_degrees,(homology_degrees[i]+1,homology_degrees[i]+2))
    end
end
for pair in consecutive_differential_degrees
    k = pair[2]
    d_k_minus_1 = LowCohomologySOS.embed.(identity, d[k-1], Ref(RG_Δ[k-1]))
    d_k = LowCohomologySOS.embed.(identity, d[k], Ref(RG_Δ[k-1]))
    @assert d_k_minus_1*d_k == [zero(RG_Δ[k-1]) for i in 1:cells_number[k-2],j in 1:cells_number[k]]
end

# The stabiliser parts which we have to add to get free modules.
# Hopefully the stabilisers' elements belong to the half bases.
stab_part_dim = Dict()
for k in homology_degrees
    stab_part_dim[k] = [
        i == j ? one(RG_Δ[k])-SLnCohomology.averaged_rep(m_stabs[k][i],RG_Δ[k]) : zero(zero(RG_Δ[k])) 
        for i in 1:cells_number[k],j in 1:cells_number[k]
    ]
end

# Compute the Laplacians.
# At the end, we embed the Laplacian into RG_star, the group ring
# with the same basis as RG but with twisted multiplciation, i.e.
# (1+g)(1+h)=1+g+h+g^(-1)h. This is needed to translate hermitian squares
# to the standard ones for the definition of the semi-definite optimization problem
# (solvers prefer standard squares to hermitian :)).
# This has no effect on the shape of the Laplacian since we just embed it.
RG_Δ_star = Dict()
for k in homology_degrees
    RG_Δ_star[k] = LowCohomologySOS.group_ring(sln, half_basis_Δ[k], star_multiplication = true)
end
Δ = Dict()
for pair in consecutive_differential_degrees
    k = pair[1]
    d_k = LowCohomologySOS.embed.(identity, d[k], Ref(RG_Δ[k]))
    d_k_plus_1 = LowCohomologySOS.embed.(identity, d[k+1], Ref(RG_Δ[k]))
    Δ[k] = d_k'*d_k+d_k_plus_1*d_k_plus_1'+stab_part_dim[k]
    Δ[k] = LowCohomologySOS.embed.(identity, Δ[k], Ref(RG_Δ_star[k]))
end
if homology_degrees[1] == differential_degrees[1]
    kx = homology_degrees[1]
    d_kx_plus_1 = LowCohomologySOS.embed.(identity, d[kx+1], Ref(RG_Δ[kx]))
    Δ[kx] = d_kx_plus_1*d_kx_plus_1'+stab_part_dim[kx]
end
if homology_degrees[end] == differential_degrees[end]-1
    kx = differential_degrees[end]
    d_kx = LowCohomologySOS.embed.(identity, d[kx], Ref(RG_Δ[kx-1]))
    stab_part_dim[kx] = [
        i == j ? one(RG_Δ[kx-1])-SLnCohomology.averaged_rep(
            m_stabs[kx][i],RG_Δ[kx-1]) : zero(zero(RG_Δ[kx-1])
            ) 
        for i in 1:cells_number[kx],j in 1:cells_number[kx]
    ]
    Δ[kx] = d_kx'*d_kx+stab_part_dim[kx]
    push!(homology_degrees,kx)
end

# A sanity check: check if the Laplacians are hermitian.
for entry in Δ
    @assert entry[2]' == entry[2]
end

# Save the Laplacians in a serialized form in a file.
laplacian_data = Dict()
laplacian_data["laplacians"] = Δ
laplacian_data["differential_degrees"] = differential_degrees
serialize(joinpath(@__DIR__, "./precomputed_laplacians/sl"*string(n)*"_laplacians.sjl"), laplacian_data)
