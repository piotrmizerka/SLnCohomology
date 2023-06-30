# For SL(n,ℤ) matrices M1,...,Mn (stored as group elements g) and their correpsonding orientation signs σi=+-1
# compute 1/n(σ1*M1+...+σn*Mn) as an element of RG (in our case RG is the group ring of SL(n,ℤ)).
function averaged_rep(elt_list, RG)
    sum_ = sum(σ*RG(g) for (g,σ) in elt_list)
    stab_order = length(elt_list)
    return sum_//stab_order
end

# Auxililiary function describing the differentials provided by Benjamin as RG elements.
# "differential_data" contains a list of pairs (coeff,RG element) and the function returns 
# the sum of coeff*RG element for all pairs. RG element is in our case the averaged representative
# computed with "averaged_rep" funcion.
# TODO: delete this as a part of common code for slns (sl3 to adopt)
function rg_elt(differential_data)
    return sum(dd[1]*dd[2] for dd in differential_data)
end

# Representing nxn arrays as honest elements of SL(n,ℤ). 
# We apply Gaussian elimination and Euclidean algorithm.
function gelt_from_matrix(M::AbstractMatrix, sln)
    inv_result, temp_matrix = one(sln), copy(M)
    for i in 1:(N-1)
        factor, temp_matrix = reduce_column_down(i,temp_matrix, sln)
        inv_result = factor*inv_result
    end
    for i in N:-1:2
        factor, temp_matrix = reduce_column_up(i,temp_matrix, slnp_order)
        inv_result = factor*inv_result
    end
    result = inv_result^(-1)
    return result
end
function reduce_column_up(col::Integer, M::AbstractMatrix, sln)
    result_factor, result_matrix = one(sln), copy(M)
    for k in (col-1):-1:1
        mult = div(M[k,col],M[col,col])
        factor, result_matrix = add_row_multiplicity(k,col,-mult,result_matrix,sln)
        result_factor = factor*result_factor
    end
    return result_factor, result_matrix
end
function reduce_column_down(col::Integer, M::AbstractMatrix, sln)
    result_factor, result_matrix = one(sln), copy(M)
    for k in (col+1):N
        factor, result_matrix = euclidean_algorithm(col,k,col,result_matrix,sln)
        result_factor = factor*result_factor
    end
    return result_factor, result_matrix
end
function euclidean_algorithm(i_row::Integer, j_row::Integer, col::Integer, M::AbstractMatrix, sln)
    result_factor, result_matrix = one(sln), copy(M)
    it = 0
    while result_matrix[i_row,col] != 1 && result_matrix[i_row,col] != -1 && it < 10
        if abs(result_matrix[i_row,col]) <= abs(result_matrix[j_row,col]) && result_matrix[i_row,col] != 0
            mult = div(result_matrix[j_row,col],result_matrix[i_row,col])
            factor, result_matrix = add_row_multiplicity(j_row,i_row,-mult,result_matrix,sln)
            result_factor = factor*result_factor
        elseif result_matrix[j_row,col] != 0
            factor, result_matrix = flip_two_rows(i_row,j_row,result_matrix,sln)
            result_factor = factor*result_factor
        end
        it += 1
    end
    if result_matrix[j_row,col] != 0 # zeroing the jth row
        mult = div(result_matrix[j_row,col],result_matrix[i_row,col])
        factor, result_matrix = add_row_multiplicity(j_row,i_row,-mult,result_matrix,sln)
        result_factor = factor*result_factor
    end
    if result_matrix[i_row,col] != 1 # we want 1 instead of -1 in the ith row
        factor, result_matrix = invert_two_rows(i_row,j_row,result_matrix,sln)
        result_factor = factor*result_factor
    end
    return result_factor, result_matrix
end
function add_row_multiplicity(i::Integer, j::Integer, mult::Integer, M::AbstractMatrix, sln)
    E(i,j) = E_(i,j,sln)
    if mult != -1 && mult != 1
        result_factor = E(i,j)^mult
    elseif mult == 1
        result_factor = E(i,j)
    else
        result_factor = E(i,j)^(-1)
    end
    result_matrix = copy(M)
    for k in 1:N
        result_matrix[i,k] = M[i,k]+mult*M[j,k]
    end
    return result_factor, result_matrix
end
function flip_two_rows(i::Integer, j::Integer, M::AbstractMatrix, sln)
    E(i,j) = E_(i,j,sln)
    result_factor, result_matrix = E(i,j)*E(j,i)^(-1)*E(i,j), copy(M)
    for k in 1:N
        result_matrix[i,k] = M[j,k]
        result_matrix[j,k] = -M[i,k]
    end
    return result_factor, result_matrix
end
function invert_two_rows(i::Integer, j::Integer, M::AbstractMatrix, sln)
    E(i,j) = E_(i,j,sln)
    result_factor, result_matrix = E(i,j)^2*E(j,i)^(-1)*E(i,j)^2*E(j,i)^(-1), copy(M)
    for k in 1:N
        result_matrix[i,k] = -M[i,k]
        result_matrix[j,k] = -M[j,k]
    end
    return result_factor, result_matrix
end
function E_(i::Integer, j::Integer, sln)
    Eij = MatrixGroups.ElementaryMatrix{N}(i,j,Int8(1))
    for s in gens(sln)
        if MatrixGroups.matrix_repr(s) == MatrixGroups.matrix_repr(Eij)
            return s
        end
    end
end