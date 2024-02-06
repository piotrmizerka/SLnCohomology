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

# TODO: this takes too much time!
function ind_H_to_G(g, π, coset_data, deg, p::Integer)
    cosets = coset_data["cosets"]
    cosets_representatives = coset_data["cosets_representatives"]
    cosets_representatives_indices = coset_data["cosets_representatives_indices"]
    total_size = deg*length(cosets_representatives)
    result = spzeros(total_size,total_size)
    for i in eachindex(cosets_representatives)
        gi = cosets_representatives[i]
        ggi = matrix_mod_p(g*gi,p)
        gj = cosets[ggi]
        j = cosets_representatives_indices[gj]
        gj_inv = inv(AbstractAlgebra.matrix(GF(p),gj))
        h = Int8.(AbstractAlgebra.lift.(gj_inv*AbstractAlgebra.matrix(GF(p),ggi)))
        set_block(result, π[h], j, i)
    end
    return result
end

function ind_rep_dict(support, π, coset_data, deg, p::Integer)
    result = Dict()
    it = 0
    for g in support
        if it%100 == 0
            @info it, "pi dict"
        end
        result[g] = ind_H_to_G(g, π, coset_data, deg, p)
        it += 1
    end
    return result
end

function laplacians_support(Δ, p::Integer)
    support = []
    for entry in Δ
        n = entry[1]
        @info n
        RG = parent(first(Δ[n]))
        for k in eachindex(Δ[n])
            for i in SparseArrays.nonzeroinds(Δ[n][k].coeffs)
                x = RG.basis[i]
                proj = projection(MatrixGroups.matrix_repr(x),p)
                push!(support,proj)
            end
        end
    end
    return collect(Set(support))
end

function matrix_mod_p(M,p::Integer)
    result = copy(M)
    for i in eachindex(result)
        result[i] %= p
    end
    return result
end

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

function permutation_matrices(matrices, slnp_dict, p)
    result = Dict()
    for M in matrices
        left_action_matrix(i) = projection(M*matrices[i],p)
        perm = PermutationGroups.Perm([slnp_dict[left_action_matrix(i)] for i in eachindex(matrices)])
        result[M] = permutation_matrix(perm)
    end
    return result
end

function permutation_matrix(perm)
    degree = maximum(perm.d)
    result = spzeros(degree,degree)
    for j in eachindex(perm.d)
        result[perm.d[j],j] = 1
    end
    return result
end

# Compute projection onto SL(N,p) from a matrix from SL(N,Z)
function projection(M, p::Integer)
    result = Int8.(zeros(size(M)[1],size(M)[2]))
    for i in 1:size(M)[1]
        for j in 1:size(M)[2]
            result[i,j] = Int8(M[i,j]%p)
            result[i,j] = (result[i,j] >= 0 ? result[i,j] : result[i,j]+p)
        end
    end
    return result
end

function read_slnp_matrices(file_path,N::Integer)
    sl_n_p_matrices = []
    file = open(file_path, "r")
    i = 0
    current_matrix = Int8.(zeros(N,N))
    for line in eachline(file)
        linex = replace(line, r"\s+" => "")
        linexx = replace(linex, r"\." => "0")
        for j in eachindex(linexx)
            current_matrix[i%N+1,j] = parse(Int8,linexx[j])
        end
        if i%N == N-1
            push!(sl_n_p_matrices,current_matrix)
            current_matrix = Int8.(zeros(N,N))
        end
        i += 1
    end
    close(file)
    return sl_n_p_matrices
end

# Permutation matrix corresponding do a group ring elt ξ given by the perm repr
function representing_matrix(ξ, π, deg::Integer, p::Integer, subgroup_index::Integer)
    n = deg*subgroup_index
    result = typeof(first(ξ.coeffs)).(spzeros(n,n))
    RG = parent(ξ)
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ(x)*π[proj]
    end
    return Matrix(result)
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

# Save all elts from SL(n,p) - brute force approach
# Keep their order in the dictionary - we use it later to create permutation matrices
function slnp(n::Integer,p::Integer)
    result = Dict()
    values = 0:(p-1)
    tuples = collect(Iterators.product(fill(values, n^2)...))
    i = 0
    matrices_tuples_order = []
    for tuple in tuples
        candidate = [tuple[(i-1)*n+j] for i in 1:n,j in 1:n]
        det_as_integer = Int(det(candidate))
        det_mod_p = (det_as_integer%p >= 0 ? det_as_integer%p : det_as_integer%p+p)
        if det_mod_p == 1
            i += 1
            result[candidate] = i
            push!(matrices_tuples_order,candidate)
        end
    end
    @assert i == slnp_order(n,p)
    return result, matrices_tuples_order
end

function slnp_order(n,p)
    result = 1
    for i in 0:(n-1)
        result *= (p^n-p^i)
    end
    result = div(result,p-1)
    return result
end

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
