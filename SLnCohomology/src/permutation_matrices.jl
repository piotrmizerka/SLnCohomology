function slnp_order(n,p)
    result = 1
    for i in 0:(n-1)
        result *= (p^n-p^i)
    end
    result = div(result,p-1)
    return result
end

# Save all elts from SL(n,p) - brute force approach
function slnp(n::Integer,p::Integer)
    result = Dict()
    values = 0:(p-1)
    tuples = collect(IterTools.product(fill(values, n^2)...))
    i = 0
    matrices_tuples_order = []
    for tuple in tuples
        candidate = [tuple[(i-1)*n+j] for i in 1:n,j in 1:n]
        det_as_integer = Int(round(det(candidate)))
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
    degree =maximum(perm.d)
    result = spzeros(degree,degree)
    for j in eachindex(perm.d)
        result[perm.d[j],j] = 1
    end
    return result
end
