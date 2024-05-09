# Coset data for a subgroup H of sl_n_p_matrices:
# (1) cosets - stores representatives of each coset for every g ∈ SL(N,p),
# (2) cosets_representatives - coset representatives storedin a list,
# (3) cosets_representatives_indices - indices of cosets'representatives.
function coset_data(H, sl_n_p_matrices, p::Integer)
    cosets = Dict()
    cosets_representatives = []
    considered_matrices = Dict()
    for g in sl_n_p_matrices
        considered_matrices[g] = false
    end
    for g in sl_n_p_matrices
        if !considered_matrices[g]
            for h in H
                gh = matrix_mod_p(g*h,p)
                cosets[gh] = g
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
        "cosets" => cosets, 
        "cosets_representatives" => cosets_representatives, 
        "cosets_representatives_indices" => cosets_representatives_indices
    )
end

# Saves, for each element of support, induced representation in a dictionary.
function ind_rep_dict(support, π, coset_data, deg::Integer, p::Integer)
    result = Dict()
    for g in support
        result[g] = ind_H_to_G(g, π, coset_data, deg, p)
    end
    return result
end

# Support of Laplacians (common for all o them).
function laplacians_support(Δ, p::Integer)
    support = []
    for entry in Δ
        n = entry[1]
        RG = parent(first(Δ[n]))
        for k in eachindex(Δ[n])
            for i in SparseArrays.nonzeroinds(Δ[n][k].coeffs)
                x = RG.basis[i]
                proj = matrix_mod_p(MatrixGroups.matrix_repr(x),p)
                push!(support,proj)
            end
        end
    end
    return collect(Set(support))
end

# Check a sufficient condition for non-existence of invariant vectors.
# Returns a subspace which is always non-invariant.
# Yields sth nontrivial for flip-perm reps only.
function no_inv_subspace(H, π)
    no_inv_subspace = Set([])
    for h in H
        no_inv_h = []
        for i in 1:size(π[h])[1]
            if π[h][i,i] == -1
                push!(no_inv_h,i)
            end
        end
        no_inv_subspace = union!(no_inv_subspace,no_inv_h)
    end
    return no_inv_subspace
end

# Maximal degree of permutations from a given list.
function permutations_degree(permutations)
    deg = 0
    for perm in permutations
        for cycle in perm
            if length(cycle) > 0
                if deg < maximum(cycle)
                    deg = maximum(cycle)
                end
            end
        end
    end
    return deg
end

# Permutation matrix corresponding to a group ring elt ξ given by the repr π
function representing_matrix(ξ, π, p::Integer)
    RG = parent(ξ)
    n = size(π[matrix_mod_p(MatrixGroups.matrix_repr(first(RG.basis)),p)])[1]
    result = typeof(first(ξ.coeffs)).(spzeros(n,n))
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = matrix_mod_p(MatrixGroups.matrix_repr(x),p)
        result += ξ(x)*π[proj]
    end
    return Matrix(result)
end

# Save a permutation in a more suitable format for further parsing.
function standarize_permutation(perm, deg)
    result = (perm[1] == [] ? [] : perm)
    considered = [false for i in 1:deg]
    for cycle in perm
        for x in cycle
            considered[x] = true
        end
    end
    missing_numbers = []
    for i in 1:deg
        if !considered[i]
            push!(missing_numbers,i)
        end
    end
    for x in missing_numbers
        push!(result,[x])
    end
    return result
end
