# For SL(n,ℤ) matrices M1,...,Mn (stored as group elements g) and their correpsonding orientation signs σi=+-1
# compute 1/n(σ1*M1+...+σn*Mn) as an element of RG (in our case RG is the group ring of SL(n,ℤ)).
function averaged_rep(elt_list, RG)
    sum_ = sum(σ*RG(g) for (g,σ) in elt_list)
    stab_order = length(elt_list)
    return sum_//stab_order
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

# Permutation matrix of a permutation.
# Used only for the regular rep of SL(3,3) which has been to moved tests.
function permutation_matrix(perm)
    degree = maximum(perm.d)
    result = spzeros(Int64,degree,degree)
    for j in eachindex(perm.d)
        result[perm.d[j],j] = 1//1
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
    for it in SparseArrays.nonzeroinds(sparse(vec(result)))
        @assert typeof(result[it]) <: Rational
    end
    return Matrix(result)
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
        temp = detx(candidate)%p
        det_mod_p = (temp >= 0 ? temp : temp+p)
        if det_mod_p == 1
            push!(result,candidate)
        end
    end
    return result
end

function slnp_order(n::Integer, p::Integer)
    result = 1
    for i in 0:(n-1)
        result *= (p^n-p^i)
    end
    result = div(result,p-1)
    return result
end
