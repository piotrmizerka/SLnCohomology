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
# The degrees of the Laplacian that the user wants to compute
user_demanded_degrees = [parse(Int64, ARGS[i]) for i in 2:length(ARGS)]

sln = MatrixGroups.SpecialLinearGroup{n}(Int8)

# The boundary and stabiliser data
cells_sln = SLnCohomology.cells_sln(n)
oriented_cells_sln = SLnCohomology.oriented_cells_dict(cells_sln)
boundaries = SLnCohomology.boundaries_dict(oriented_cells_sln)
stabilisers = SLnCohomology.stabilisers_dict(oriented_cells_sln)

@info "Voronoi tesselation computed"

cell_degrees = sort([first(x) for x in boundaries]) # the degrees where the Voronoi complex has cells intersecting the interior
min_degree = cell_degrees[1]
max_degree = cell_degrees[end]

# Remove demanded degrees that don't have non-trivial cells
reasonable_demanded_degrees = []
for degree in user_demanded_degrees
    if degree in cell_degrees
        push!(reasonable_demanded_degrees, degree)
    end
end

relevant_degrees_differential = [] # to compute the Laplacian Delta_n, we need the differentials partial_n and partial_{n+1}
for degree in cell_degrees
    if (degree in reasonable_demanded_degrees)||(degree-1 in reasonable_demanded_degrees)
        push!(relevant_degrees_differential, degree)
    end
end
relevant_degrees_group_ring = [] # we need the group rings in one degree higher, as the images of the adjoints
for degree in cell_degrees
    if (degree in reasonable_demanded_degrees)||(degree-1 in reasonable_demanded_degrees)||(degree+1 in reasonable_demanded_degrees)
        push!(relevant_degrees_group_ring, degree)
    end
end

# Compute the supports (i.e. half_bases) for the group rings to compute the Laplacians
# - we just add to half_basis the coset elements appearing in the differentials.
d_union = Dict(k=>[one(sln)] for k in cell_degrees)
for k in relevant_degrees_group_ring
    if !(k == min_degree) # in min_degree, all boundaries are trivial, so nothing needs to be added
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
end
half_basis_Δ = Dict()
for k in relevant_degrees_group_ring
    if k == max_degree
        half_basis_Δ[k] = d_union[k]
        half_basis_Δ[k] = unique([half_basis_Δ[k];inv.(half_basis_Δ[k])])
    else
        half_basis_Δ[k] = union(d_union[k], d_union[k+1])
        half_basis_Δ[k] = unique([half_basis_Δ[k];inv.(half_basis_Δ[k])])
    end
end

# Compute the group rings (with standard multiplication by convolution: (1+g)(1+h)=1+g+h+gh)
RG_Δ = Dict()
for k in relevant_degrees_group_ring
    RG_Δ[k] = LowCohomologySOS.group_ring(sln, half_basis_Δ[k])
    @info "Computed group ring support in degree $k"
end

# Compute the differentials as matrices over RGs.
# We store them as matrices over Rationals to perform exact operations.
cells_number = Dict(k=>length(cells_sln[k]) for k in cell_degrees) # number of cell orbits in degree k
cells_number[min_degree-1] = 0
dx = Dict([(k,i,j) => 0//1*zero(RG_Δ[k-1]) for k in relevant_degrees_differential for i in 1:cells_number[k-1] for j in 1:cells_number[k]]) # This is empty for min_degree
for k in relevant_degrees_differential
    for j in eachindex(boundaries[k]) # each j corresponds to a k-cell sigma
        for tauprime in boundaries[k][j] # each tau' is a facet of sigma
            coset_orient = [(SLnCohomology.gelt_from_matrix(M,sln),eta) for (M,eta) in tauprime["orbit_coset_with_orientation"]]
            # coset_orient is G(tau,tau'); gelt_from_matrix(M,sln) is an element g sending tau to tau', eta is eta(tau,tau',g)
            summand_in_RG = tauprime["sign"]*SLnCohomology.averaged_rep(coset_orient,RG_Δ[k-1]) 
            # summand_in_RG is x_{tau'}; tauprime["sign"] is eta(tau',sigma,1)
            dx[k,tauprime["orbit_standard_cell"],j] += summand_in_RG 
            # dx[k,tauprime["orbit_standard_cell"],j] is partial_{sigma,tau}; tauprime["orbit_standard_cell"] is tau
        end
    end 
end
d = Dict()
for k in relevant_degrees_differential
    d[k] = [dx[k,i,j] for i in 1:cells_number[k-1], j in 1:cells_number[k]]
end

# Sanity checks for vanishing compositions of differentials.
consecutive_relevant_degrees = []
for i in relevant_degrees_differential
    if i+1 in reasonable_demanded_degrees
        push!(consecutive_relevant_degrees,(i,i+1))
    end
end
for pair in consecutive_relevant_degrees
    k = pair[2]
    d_k_minus_1 = LowCohomologySOS.embed.(identity, d[k-1], Ref(RG_Δ[k-1]))
    d_k = LowCohomologySOS.embed.(identity, d[k], Ref(RG_Δ[k-1]))
    @assert d_k_minus_1*d_k == [zero(RG_Δ[k-1]) for i in 1:cells_number[k-2],j in 1:cells_number[k]]
end

# The stabiliser parts which we have to add to get free modules.
# Convert the stabilisers to signed subgroups of SL(n,ℤ)
stab_part_dim = Dict()
for k in reasonable_demanded_degrees
    m_stabs = [
        [(SLnCohomology.gelt_from_matrix(M,sln),eta) for (M,eta) in stab] for stab in stabilisers[k]
    ]
    stab_part_dim[k] = [
        i == j ? one(RG_Δ[k])-SLnCohomology.averaged_rep(m_stabs[i],RG_Δ[k]) : zero(zero(RG_Δ[k])) 
        for i in 1:cells_number[k],j in 1:cells_number[k]
    ]
    # The above also verifies that the stabilisers' elements belong to the half bases.
end

# LEt's add a test saying that x^*x = x for x =  one(RG_Δ[k])-SLnCohomology.averaged_rep(m_stabs[i],RG_Δ[k]) : zero(zero(RG_Δ[k]))

# Compute the Laplacians. (Currently in `meta code' -- to be corrected, as in line 141)
Δ = Dict()
for k in reasonable_demanded_degrees
    if k == min_degree
        d_k_plus_1 = LowCohomologySOS.embed.(identity, d[k+1], Ref(RG_Δ[k]))
        Δ[k] = d_k_plus_1*copy(d_k_plus_1')*[ID - stab_part_dim[k]] +stab_part_dim[k]
    elseif k == max_degree
        d_k = LowCohomologySOS.embed.(identity, d[k], Ref(RG_Δ[k]))
        Δ[k] = copy(d_k')*d_k+stab_part_dim[k]
    else
        d_k = LowCohomologySOS.embed.(identity, d[k], Ref(RG_Δ[k]))
        d_k_plus_1 = LowCohomologySOS.embed.(identity, d[k+1], Ref(RG_Δ[k]))
        Δ[k] = copy(d_k')*d_k+d_k_plus_1*copy(d_k_plus_1')+stab_part_dim[k]
    end
end

# A sanity check: check if the Laplacians are hermitian.
for entry in Δ
    @assert entry[2]' == entry[2]
end

# Save the Laplacians in a serialized form in a file.
serialize(joinpath(@__DIR__, "precomputed_laplacians/sl"*string(n)*"_laplacians.sjl"), Δ)