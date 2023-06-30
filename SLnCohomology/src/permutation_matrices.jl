function slnp_order(n,p)
    result = 1
    for i in 0:(n-1)
        result *= (p^n-p^i)
    end
    result = div(result,p-1)
    return result
end

# Save all elts from SL(n,p) - brute force approach
# Keep their order in the dictionary - we use it later to create permutation matrices
function slnp(n::Integer,p::Integer)
    result = Dict()
    values = 0:(p-1)
    tuples = collect(IterTools.product(fill(values, n^2)...))
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

function permutation_matrix(perm)
    degree = maximum(perm.d)
    result = spzeros(degree,degree)
    for j in eachindex(perm.d)
        result[perm.d[j],j] = 1
    end
    return result
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

# Permutation matrix corresponding do a group ring elt ξ given by the perm repr
function representing_matrix(ξ,p::Integer,perm_mats)
    coeff_type = typeof(first(ξ.coeffs))
    deg = size(first(perm_mats)[2])[1]
    result = sparse(coeff_type.(zeros(deg,deg)))
    RG = parent(ξ)
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ.coeffs[i]*perm_mats[proj]
    end
    return result
end
