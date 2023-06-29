using LinearAlgebra
using Combinatorics
using Serialization

function signed_matrices(dimension)
    # Returns all diagonal square "sign matrices" in dimension
    signs_list = []
    for set in powerset(1:dimension)
        matrix = zeros(Int8,dimension,dimension)
        for i in 1:dimension
            if i in set
                matrix[(i-1)*dimension+(i)] = -1
            else
                matrix[(i-1)*dimension+(i)] = 1
            end
        end
        push!(signs_list,matrix)
    end    
    return signs_list
end

function bases_in_integer_lattices(matrix)
    #= matrix an integer matrix whose columns form a basis of the lattice
    returns all sets of bases, with all possible sign combinations 
    =#
    bases = []
    columns = collect(eachcol(matrix)) # set of all columns
    dimension = length(columns[1])
    sign_changes = signed_matrices(dimension)
    for candidate in Combinatorics.permutations(columns,dimension) # all subsets of columns of with dimension-many elements
        matrix_candidate = transpose(reduce(vcat,transpose.(candidate)))
        if (round(det(matrix_candidate)) == 1) || (round(det(matrix_candidate)) == -1)
            #if that's the case, candidate is a basis of Z^dimension
            for sign in sign_changes
                # each basis element can be multiplied by +-1. Add all combinations
                push!(bases, matrix_candidate*sign)
            end
        end
    end
    return bases
end

function relative_orientation(matrix1,matrix2)
    #= matrix1, matrix2 need to be integer matrices such that the set of lines spanned by the columns of 
        matrix1 is the same as that of matrix2
        computes the orientation of the simplex spanned by the columns
    =# 
    permutation = Int64[]
    for vector in eachcol(matrix2)
        vertex_index_in_matrix1 = 1
        vector_is_column_of_matrix1 = false
        for vector_standard in eachcol(matrix1)
            if vector == vector_standard || -vector == vector_standard
                push!(permutation,vertex_index_in_matrix1)
                vector_is_column_of_matrix1 = true
            end
            vertex_index_in_matrix1 += 1
        end
        vector_is_column_of_matrix1 || error("matrix1 isn't a permutation of matrix2.")
    end
    sign = Combinatorics.levicivita(permutation)
    return sign
end

function quadratic_form(matrix)
    form = matrix*transpose(matrix)
    return form
end

function stabiliser_coset_with_orientation(matrix1, matrix2)
    #= 
    matrix1, matrix2 integer matrix whose columns form a basis of the lattice.
    Returns all (g,orientation) such that g\in SL_dimension(Z) sends the quadratic form associated to matrix1 to 
    the one associated to matrix 2 and orientation is the relative orientation of the corresponding simplices
    =#
    elements = []
    dimension = size(matrix1)[1]
    quadratic_form1 = quadratic_form(matrix1) # the quadratic form associated to the lattice
    quadratic_form2 = quadratic_form(matrix2) # the quadratic form associated to the lattice
    bases_in_image = bases_in_integer_lattices(matrix2) # bases in the image
    basis_in_matrix_1 = bases_in_integer_lattices(matrix1)[1] # a fixed basis in the domain
    basis_in_matrix_1_inverse = Int.(basis_in_matrix_1^(-1))  #inverse as integer matrix
    for basis in bases_in_image
        g_transpose = basis*basis_in_matrix_1_inverse # this sends the fixed basis_in_matrix_1 to basis in matrix2
        if round(det(g_transpose)) == 1
            # then it lies in SL
            g = transpose(g_transpose)
            if (g_transpose*quadratic_form1*g == quadratic_form2) # then g sends the first form to the second one
                orientation = relative_orientation(g_transpose*matrix1,matrix2)
                # this is the relative orientation of g.first cell compared to second cell
                push!(elements, (g,orientation))
            end
        end
    end
    return elements
end

function same_orbit(matrix1,matrix2)
    #= Checks whether matrix1 lies in the SL_n orbit as matrix2
    =#
    if length(stabiliser_coset_with_orientation(matrix1,matrix2)) > 0
        return true
    else
        return false
    end
end

function orbit_in_list(matrix,list)
    #= Checks whether matrix lies in the SL_n orbit of one of the elements in list
    =#
    for matrix2 in list
        if same_orbit(matrix,matrix2)
            return true
        end
    end
    return false
end

function cell_in_interior(matrix)
    #= check whether the matrix determines a cell in the interior of the symmetric space
    =#
    if length(bases_in_integer_lattices(matrix))>0
        return true
    else
        return false
    end
end

function boundaries_in_group_ring_with_orientation(matrix,cell_dictionary)
    #= cell dictionary should have the cells of the complex
    This only works if the cells are simplicial!
    =#
    boundary_cells = []
    cell_dimension = size(matrix)[2] - 1 # only true if the cell is simplicial
    for j in 1:size(matrix)[2]
        sign = (-1)^(j-1) # the sign of the differential
        boundary = matrix[begin:end,begin:end .!=j] # exclude the j-th column
        if cell_in_interior(boundary)
            boundary_information = Dict()
            boundary_information["sign"] = sign
            for cell_index in 1:length(cell_dictionary[cell_dimension-1])
                # check in the list of codimension 1 cells which orbit we have
                standard_cell = cell_dictionary[cell_dimension-1][cell_index]
                if same_orbit(standard_cell,boundary)
                    boundary_information["orbit_standard_cell"] = cell_index
                    boundary_information["orbit_coset_with_orientation"] = stabiliser_coset_with_orientation(standard_cell,boundary)
                end
            end
            push!(boundary_cells,boundary_information)
        end
    end
    return boundary_cells
end