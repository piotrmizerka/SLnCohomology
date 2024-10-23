function r_pinv(A)
    #= Computes the rational pseudo inverse of an integral matrix.
    Probably slower than the usual pinv function, but uses exact computations.
    =#
    if !(typeof(A) == Matrix{Int64} || typeof(A) == Matrix{Rational})
        error("The input needs to be a rational matrix")
    end
    r_A = Rational.(A) # make rational matrix to avoid promotion to float
    return (r_A'*r_A)\r_A'
end

function boundaries_dict(cell_dictionary)
    #= computes the data of the differential associated to cell_dictionary
    =#
    boundaries_sln = Dict()
    for (dimension, cell_list) in cell_dictionary
        boundaries_this_dimension = []
        for (cell,basis) in cell_list
            push!(boundaries_this_dimension,boundaries_in_group_ring_with_orientation(cell,basis,dimension,cell_dictionary))
        end
        boundaries_sln[dimension] = boundaries_this_dimension
    end
    return boundaries_sln
end

function boundaries_in_group_ring_with_orientation(cell,basis,cell_dimension,cell_dictionary)
    #= computes the facets of cell; 
        for each facet, determines whether it lies in the orbit of a cell from cell_dictionary;
        determines all elements from SL_n(Z) sending the facet to this "standard cell".
        cell_dictionary should have the standard cells of the complex, together with a choice of 
        orientation of each such standard cell
    =#
    boundary_cells = []
    lowest_dim = minimum(keys(cell_dictionary))
    # Create a list of the numbers of vertices of cells in codimension 1
    if cell_dimension == lowest_dim
        # only trivial boundaries for cells of lowest dimension
        return boundary_cells
    end
    # determine all possible sizes of vertex sets of codim-1 cells
    sizes_in_codimension_1 = Set()
    for (codim1_cell,codim1_basis) in cell_dictionary[cell_dimension-1]
        push!(sizes_in_codimension_1,size(codim1_cell)[2])
    end
    
    vertices = collect(eachcol(cell)) # set of vertices of cell
    for vertex_number in sizes_in_codimension_1
        for potential_face in Combinatorics.combinations(vertices,vertex_number)
            #look at all subsets of the vertex set that could form a cell in codimension 1
            potential_face_matrix = transpose(reduce(vcat,transpose.(potential_face)))
            # find a basis for that face
            forms_face = []
            for minimal_vector in eachcol(potential_face_matrix)
                push!(forms_face, vec(quadratic_form(minimal_vector)))
            end
            basis_face = extract_basis(forms_face)

            # Now find out whether it is in the orbit of a standard cell in codimension 1
            for cell_index in 1:length(cell_dictionary[cell_dimension-1])
                # check in the list of codimension 1 cells which orbit we have
                (standard_cell,standard_basis) = cell_dictionary[cell_dimension-1][cell_index]
                if size(standard_cell)[2] == size(potential_face_matrix)[2]
                    # then they have the same number of vertices, now check whether they lie in the same orbit
                    # this is done by checking whether their barycenters lie in the same orbit
                    stab_coset = stabiliser_coset_with_orientation((standard_cell,standard_basis),(potential_face_matrix,basis_face))
                    if length(stab_coset) != 0
                        for form in eachcol(basis)
                            if !(form in collect(eachcol(basis_face)))
                                global extended_basis = hcat(basis_face,form)
                                break
                            end
                        end
                        sign = relative_orientation_bases(basis, extended_basis) #this is eta(tau',sigma,1)
                        boundary_information = Dict()
                        boundary_information["sign"] = sign
                        boundary_information["orbit_standard_cell"] = cell_index
                        boundary_information["orbit_coset_with_orientation"] = stab_coset
                        push!(boundary_cells,boundary_information)
                        break
                    end
                end
            end
        end
    end
    return boundary_cells
end

function extract_basis(list_of_vectors)
    #= Picks a linearly independent subset of list_of_vectors (required to be integral or rational)
    =#
    partial_basis = hcat(list_of_vectors[1])
    for vector in list_of_vectors # slightly wasteful to check first element again
        potential_partial_basis = hcat(partial_basis,vector)
        if rankx(potential_partial_basis) > rankx(partial_basis)
            partial_basis = potential_partial_basis
        end
    end
    return partial_basis
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

function oriented_cells_dict(unoriented_cell_dictionary)
    #= adds an orientation to every cell in unoriented_cell_dictionary
    =#
    cell_dictionary = Dict()
    for dimension in keys(unoriented_cell_dictionary)
        cell_dictionary[dimension] = []
        for cell in unoriented_cell_dictionary[dimension]
            forms_cell = []
            for minimal_vector in eachcol(cell)
                push!(forms_cell, vec(quadratic_form(minimal_vector)))
            end
            basis = extract_basis(forms_cell)
            push!(cell_dictionary[dimension],(cell,basis))
        end
    end
    return cell_dictionary   
end

function quadratic_form(matrix)
    #= Determines the quadratic form given by the sum of the rank-1 forms associated to each
        column of matrix.
    =#
    return matrix*transpose(matrix)
end

function relative_orientation_bases(basis1,basis2)
    #= Assumes that both basis1, basis2 are ordered basis of the same vector space.
        Computes the relative orientation
    =#
    base_change = r_pinv(basis1)*basis2 # matrix sending basis1 to basis2 in the corresponding subspace 
                                        # (expressed in basis1 of that subspace)
    if detx(base_change)>0
        return 1
    end
    return -1
end

function same_orbit(matrix1,matrix2)
    #= Checks whether matrix1 lies in the SL_n orbit as matrix2
    =#
    if length(stabiliser_coset_SL(quadratic_form(matrix1),quadratic_form(matrix2))) > 0
        return true
    end
    return false
end

function stabiliser_coset_SL(form1, form2)
    #=  Computes all elements in SL_n that send form1 to form2
    =#
    stabiliser = []
    # could add a determinant check here to speed up computation if needed
    max_norm = maximum(diag(form2)) # max diagonal entry of form2 = max form1-norm of image of a vector under g
    short_vectors1 = shortestVectors(form1,max_norm)
    for g in pleskenSouvignier_two_matrices(form1,form2,short_vectors1)
        if detx(g) == 1
            push!(stabiliser, g)
        end
    end
    return stabiliser
end

function stabiliser_coset_with_orientation((matrix1,oriented_basis1), (matrix2,oriented_basis2))
    #= 
    matrix1, matrix2 integer matrix whose columns form a basis of the lattice.
    oriented_basis1 an orientation of the form associated to matrix1, oriented_basis2 same for matrix2
    Returns all (g,orientation) such that g\in SL_n(Z) sends the quadratic form associated to matrix1 to 
    the one associated to matrix 2 and orientation is the relative orientation of the corresponding oriented cells
    =#
    n = size(matrix1)[1]
    elements = []
    for g in stabiliser_coset_SL(quadratic_form(matrix1), quadratic_form(matrix2))
        # these are the elements in the stabiliser
        # now determine the orientation
        # first compute what g does with oriented_basis1
        g_basis1 = Array{Int64}(undef, n*n, 0) # an n^2x0 matrix, to be filled with the basis
        for column in eachcol(oriented_basis1)
            g_basis1 = hcat(g_basis1,vec(transpose(g)*reshape(column,(n,n))*g))
        end
        orientation = relative_orientation_bases(g_basis1,oriented_basis2)
        # this is eta(cell1, cell2, g), the relative orientation of cell1.g compared to cell2
        push!(elements, (g,orientation))
    end
    return elements
end

function stabilisers_dict(cell_dictionary)
    #= Computes the stabilisers of the cells in cell_dictionary
    =#
    stabilisers = Dict()
    for dimension in keys(cell_dictionary)
        println("Dimension $dimension")
        cell_list = cell_dictionary[dimension]
        stabilisers_this_dimension = []
        cell_count = 0
        for (cell,basis) in cell_list
            cell_count += 1
            println("Cell number $cell_count")
            cell_stabiliser = []
            for (g,orientation) in stabiliser_coset_with_orientation((cell,basis), (cell,basis))
                push!(cell_stabiliser,(g,orientation))
            end
            push!(stabilisers_this_dimension,cell_stabiliser)
        end
        stabilisers[dimension] = stabilisers_this_dimension
    end
    return stabilisers
end
