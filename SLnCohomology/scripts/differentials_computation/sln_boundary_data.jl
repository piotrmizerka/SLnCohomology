# Load the boundary and stabiliser data
cells_sln = SLnCohomology.cells_sln(N) # N has to be defined before including this file
oriented_cells_sln = SLnCohomology.oriented_cells_dict(cells_sln)
boundaries = SLnCohomology.boundaries_dict(oriented_cells_sln)
stabilisers = SLnCohomology.stabilisers_dict(oriented_cells_sln)

# Extract homology degrees
differential_degrees = sort([first(x) for x in boundaries])

# Compute the stabilisers as signed subgroups of SL(n,ℤ)
m_stabs = Dict()
sln = MatrixGroups.SpecialLinearGroup{N}(Int8)
for n in differential_degrees
    m_stabs[n] = [
        [(SLnCohomology.gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in stabilisers[n]
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

# Compute the differential as matrices over RGs
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