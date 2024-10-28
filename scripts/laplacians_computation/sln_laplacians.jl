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
n = parse(Int64, ARGS[1]) # add the lsit of Laplacian degreess we want to compute as a parameter
demanded_degrees = [parse(Int64, ARGS[i]) for i in 2:length(ARGS)]

# The boundary and stabiliser data
cells_sln = SLnCohomology.cells_sln(n)
oriented_cells_sln = SLnCohomology.oriented_cells_dict(cells_sln)
boundaries = SLnCohomology.boundaries_dict(oriented_cells_sln)
stabilisers = SLnCohomology.stabilisers_dict(oriented_cells_sln)

# Extract homology degrees
differential_degrees = sort([first(x) for x in boundaries])
relevant_degrees = []
for degree in differential_degrees
    if (degree in demanded_degrees)||(degree-1 in demanded_degrees)
        push!(relevant_degrees, degree)
    end
end

# Compute the stabilisers as signed subgroups of SL(n,ℤ)
m_stabs = Dict()
sln = MatrixGroups.SpecialLinearGroup{n}(Int8)
for k in relevant_degrees
    m_stabs[k] = [
        [(SLnCohomology.gelt_from_matrix(M,sln),eta) for (M,eta) in stab] for stab in stabilisers[k]
    ]
end

# Compute the supports (i.e. half_bases) for the group rings to compute the Laplacians
# - we just add to half_basis the coset elements appearing in the differentials.
d_union = Dict(k=>[one(sln)] for k in relevant_degrees[2:end])
for k in relevant_degrees[2:end]
    for i in eachindex(boundaries[k])
        for j in eachindex(boundaries[k][i])
            coset = boundaries[k][i][j]["orbit_coset_with_orientation"]
            d_union[k] = union(
                d_union[k], 
                [SLnCohomology.gelt_from_matrix(M,sln) for (M,eta) in coset]
            )
        end
    end
end
half_basis_Δ = Dict()
min_degree = relevant_degrees[1]
half_basis_Δ[min_degree] = d_union[min_degree+1]
half_basis_Δ[min_degree] = unique(
    [half_basis_Δ[min_degree];inv.(half_basis_Δ[min_degree])]
)
for k in relevant_degrees[3:end]
    half_basis_Δ[k-1] = union(d_union[k-1], d_union[k])
    half_basis_Δ[k-1] = unique([half_basis_Δ[k-1];inv.(half_basis_Δ[k-1])])
end

# Compute the group rings (with standard multiplciation by convolution: (1+g)(1+h)=1+g+h+gh)
homology_degrees = []
for k in relevant_degrees[1:end-1]
    push!(homology_degrees,k)
end
RG_Δ = Dict()
for k in homology_degrees
    RG_Δ[k] = LowCohomologySOS.group_ring(sln, half_basis_Δ[k])
end

# Compute the differentials as matrices over RGs.
# We store them as matrices over Rationals to perform exact operations.
cells_number = Dict(k=>length(stabilisers[k]) for k in relevant_degrees)
dx = Dict([(k+1,i,j) => 0//1*zero(RG_Δ[k]) for k in homology_degrees for i in 1:cells_number[k] for j in 1:cells_number[k+1]])
for k in homology_degrees
    for j in eachindex(boundaries[k+1]) # each j corresponds to a (k+1)-cell sigma
        for tauprime in boundaries[k+1][j] # each tau' is a facet of sigma
            coset_orient = [(SLnCohomology.gelt_from_matrix(M,sln),eta) for (M,eta) in tauprime["orbit_coset_with_orientation"]]
            # coset_orient is G(tau,tau'); gelt_from_matrix(M,sln) is an element g sending tau to tau', eta is eta(tau,tau',g)
            summand_in_RG = tauprime["sign"]*SLnCohomology.averaged_rep(coset_orient,RG_Δ[k]) 
            # summand_in_RG is x_{tau'}; tauprime["sign"] is eta(tau',sigma,1)
            dx[k+1,tauprime["orbit_standard_cell"],j] += summand_in_RG 
            # dx[k+1,tauprime["orbit_standard_cell"],j] is partial_{sigma,tau}; tauprime["orbit_standard_cell"] is tau
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
Δ = Dict()
for pair in consecutive_differential_degrees
    k = pair[1]
    d_k = LowCohomologySOS.embed.(identity, d[k], Ref(RG_Δ[k]))
    d_k_plus_1 = LowCohomologySOS.embed.(identity, d[k+1], Ref(RG_Δ[k]))
    Δ[k] = copy(d_k')*d_k+d_k_plus_1*copy(d_k_plus_1')+stab_part_dim[k]
end
if homology_degrees[1] == relevant_degrees[1]
    kx = homology_degrees[1]
    d_kx_plus_1 = LowCohomologySOS.embed.(identity, d[kx+1], Ref(RG_Δ[kx]))
    Δ[kx] = d_kx_plus_1*copy(d_kx_plus_1')+stab_part_dim[kx]
end
if homology_degrees[end] == relevant_degrees[end]-1
    kx = relevant_degrees[end]
    d_kx = LowCohomologySOS.embed.(identity, d[kx], Ref(RG_Δ[kx-1]))
    stab_part_dim[kx] = [
        i == j ? one(RG_Δ[kx-1])-SLnCohomology.averaged_rep(
            m_stabs[kx][i],RG_Δ[kx-1]) : zero(zero(RG_Δ[kx-1])
            ) 
        for i in 1:cells_number[kx],j in 1:cells_number[kx]
    ]
    Δ[kx] = copy(d_kx')*d_kx+stab_part_dim[kx]
    push!(homology_degrees,kx)
end

# A sanity check: check if the Laplacians are hermitian.
for entry in Δ
    @assert entry[2]' == entry[2]
end

# Save the Laplacians in a serialized form in a file.
laplacian_data = Dict()
laplacian_data["laplacians"] = Δ
laplacian_data["differential_degrees"] = relevant_degrees
serialize(joinpath(@__DIR__, "./precomputed_laplacians/sl"*string(n)*"_laplacians.sjl"), laplacian_data)
