# Coset data for a subgroup H of G (G is assumed to be a subgroup of SL(N,p)):
# (1) elt_coset_labels - stores representatives of each coset for every g ∈ G,
# (2) cosets_representatives - coset representatives storedin a list,
# (3) cosets_representatives_indices - indices of cosets'representatives.
function coset_data(H, G, p::Integer)
    elt_coset_labels = Dict()
    cosets_representatives = []
    considered_matrices = Dict()
    for g in G
        considered_matrices[g] = false
    end
    for g in G
        if !considered_matrices[g]
            for h in H
                gh = matrix_mod_p(g*h,p)
                elt_coset_labels[gh] = g
                considered_matrices[gh] = true
            end
            push!(cosets_representatives,g)
        end
    end
    cosets_representatives_indices = Dict()
    for i in eachindex(cosets_representatives)
        cosets_representatives_indices[cosets_representatives[i]] = i
    end
    return Dict(
        "elt_coset_labels" => elt_coset_labels, 
        "cosets_representatives" => cosets_representatives, 
        "cosets_representatives_indices" => cosets_representatives_indices
    )
end

# Representation matrix of g∈G induced from the rep π of H≤G.
# Requires also providing coset_data, degree of π, and modulus p.
function ind_H_to_G(g, π, coset_data, deg::Integer, p::Integer)
    elt_coset_labels = coset_data["elt_coset_labels"]
    cosets_representatives = coset_data["cosets_representatives"]
    cosets_representatives_indices = coset_data["cosets_representatives_indices"]
    total_size = deg*length(cosets_representatives)
    result = spzeros(Int64,total_size,total_size)
    for i in eachindex(cosets_representatives)
        gi = cosets_representatives[i]
        ggi = matrix_mod_p(g*gi,p)
        gj = elt_coset_labels[ggi]
        j = cosets_representatives_indices[gj]
        gj_inv = inv(AbstractAlgebra.matrix(GF(p),gj))
        h = Int8.(AbstractAlgebra.lift.(gj_inv*AbstractAlgebra.matrix(GF(p),ggi)))
        set_block(result, π[h], j, i)
    end
    return result
end

# Saves, for each element of support, induced representation in a dictionary.
function ind_rep_dict(support, π, coset_data, deg::Integer, p::Integer)
    result = Dict()
    for g in support
        result[g] = ind_H_to_G(g, π, coset_data, deg, p)
    end
    return result
end

# Fills (i,j)th square block of a square matrix M with a square matrix B  
function set_block(M, B, i::Integer, j::Integer)
    deg = size(B)[1]
    for it in SparseArrays.nonzeroinds(sparse(vec(B)))
        reminder, ratio = it%deg, div(it,deg)
        k = (reminder == 0) ? deg : reminder
        l = (it%deg == 0) ? ratio : (ratio+1)
        M[(i-1)*deg+k,(j-1)*deg+l] = B[k,l]
    end
end