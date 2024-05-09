# For SL(n,ℤ) matrices M1,...,Mn (stored as group elements g) and their correpsonding orientation signs σi=+-1
# compute 1/n(σ1*M1+...+σn*Mn) as an element of RG (in our case RG is the group ring of SL(n,ℤ)).
function averaged_rep(elt_list, RG)
    sum_ = sum(σ*RG(g) for (g,σ) in elt_list)
    stab_order = length(elt_list)
    return sum_//stab_order
end

# Compute the determinant of a integer matrix.
# We don't use the built-in det function directly as it converts
# determinants of integer matrices to floats by default 
# and we want to work with exact numbers. Instead, we use the
# built-in det after casting to rationals to ensure exact operations.
function determinant(M)
    return Int8(det([M[i,j]//1 for i in 1:size(M)[1],j in 1:size(M)[2]]))
end

# Representation matrix of g∈G induced from the rep π of H≤G.
# Requires also providing coset_data, degree of π, and modulus p.
function ind_H_to_G(g, π, coset_data, deg::Integer, p::Integer)
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

# Compute projection onto SL(N,p) from a matrix from SL(N,Z)
function matrix_mod_p(M,p::Integer)
    result = Int8.(zeros(size(M)[1],size(M)[2]))
    for i in 1:size(M)[1]
        for j in 1:size(M)[2]
            result[i,j] = Int8(M[i,j]%p)
            result[i,j] = (result[i,j] >= 0 ? result[i,j] : result[i,j]+p)
        end
    end
    return result
end

# Permutation matrices of permutations correpsonding the regular rep of matrices from SL(N,p).
# Used only for the regular rep of SL(3,3) which has been to moved tests.
function regular_rep(matrices, slnp_dict, p::Integer)
    result = Dict()
    for M in matrices
        left_action_matrix(i) = matrix_mod_p(M*matrices[i],p)
        perm = PermutationGroups.Perm([slnp_dict[left_action_matrix(i)] for i in eachindex(matrices)])
        result[M] = permutation_matrix(perm)
    end
    return result
end

# Permutation matrix of a permutation.
# Used only for the regular rep of SL(3,3) which has been to moved tests.
function permutation_matrix(perm)
    degree = maximum(perm.d)
    result = spzeros(degree,degree)
    for j in eachindex(perm.d)
        result[perm.d[j],j] = 1
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

# Store all elements of SL(n,p) as matrices over integers.
function sl_n_p(
    N::Integer,
    p::Integer
)
    result = []
    all_tuples = collect(Iterators.product(fill(0:(p-1), N^2)...))
    for tuple in all_tuples
        candidate = ones(Int8,N,N)
        for it in eachindex(tuple)
            i, j = div(it-1,N)+1, (it-1)%N+1
            candidate[i,j] = tuple[it]
        end
        temp = determinant(candidate)%p
        det_mod_p = (temp >= 0 ? temp : temp+p)
        if det_mod_p == 1
            push!(result,candidate)
        end
    end
    return result
end

# Save all elts from SL(n,p) - brute force approach
# Keep their order in the dictionary - we use it later to create permutation matrices
# Used only for the regular rep of SL(3,3) which has been to moved tests.
function slnp(n::Integer,p::Integer)
    result = Dict()
    values = 0:(p-1)
    tuples = collect(Iterators.product(fill(values, n^2)...))
    i = 0
    matrices_tuples_order = []
    for tuple in tuples
        candidate = [tuple[(i-1)*n+j] for i in 1:n,j in 1:n]
        det_as_integer = round(Int,det(candidate))
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

function slnp_order(n::Integer, p::Integer)
    result = 1
    for i in 0:(n-1)
        result *= (p^n-p^i)
    end
    result = div(result,p-1)
    return result
end
