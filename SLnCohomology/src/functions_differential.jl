using LinearAlgebra
using Combinatorics
using Serialization
include("Plesken-Souvignier.jl");

function quadratic_form(matrix)
    form = matrix*transpose(matrix)
    return form
end

function relative_orientation_bases(basis1,basis2)
    #= Assume that both basis1, basis2 are ordered basis of the same vector space. 
        Computes the relative orientation
    =#
    base_change = pinv(basis1)*basis2 # matrix sending basis1 to basis2 in the corresponding subspace
    # has no check for det=0 yet, so assumes correct input
    if det(base_change)>0
        return 1
    else
        return -1
    end
end

function extract_basis(list_of_vectors)
    #= Picks a linearly independent subset of list_of_vectors
    =#
    partial_basis = hcat(list_of_vectors[1])
    for vector in list_of_vectors # bit wasteful to check first element again
        potential_partial_basis = hcat(partial_basis,vector)
        if rank(potential_partial_basis)> rank(partial_basis)
            partial_basis = potential_partial_basis
        end
    end
    return partial_basis
end

function stabiliser_coset_SL(form1, form2)
    # Compute all elements in SL_n that send form1 to form2
    stabiliser = []
    
    if round(det(form1)) !== round(det(form2)) 
        # if they don't have the same determinant, they can't lie in the same orbit
        # just to speed up a bit
        return stabiliser
    end
    
    max_norm = maximum(diag(form2)) # max diagonal entry of form2 = max form1-norm of image of a vector under g
    short_vectors1 = shortestVectors(form1,max_norm)
    for g in pleskenSouvignier_two_matrices(form1,form2,short_vectors1)
        if round(det(g)) == 1
            push!(stabiliser, g)
        end
    end
    return stabiliser
end

function stabiliser_coset_with_orientation((matrix1,oriented_basis1), (matrix2,oriented_basis2))
    #= 
    matrix1, matrix2 integer matrix whose columns form a basis of the lattice.
    oriented_basis1 an orientation of the form associated to matrix1, oriented_basis2 same for matrix2
    Returns all (g,orientation) such that g\in SL_dimension(Z) sends the quadratic form associated to matrix1 to 
    the one associated to matrix 2 and orientation is the relative orientation of the corresponding oriented cells
    =#
    n = size(matrix1)[1]
    elements = []
    for g in stabiliser_coset_SL(quadratic_form(matrix1), quadratic_form(matrix2))
        # these are the elements in the stabiliser
        # now need to determine the orientation
        # first compute what g does with oriented_basis1
        g_basis1 = Array{Int64}(undef, n*n, 0) # an n^2x0 matrix, to be filled with the basis
        for column in eachcol(oriented_basis1)
            g_basis1 = hcat(g_basis1,vec(transpose(g)*reshape(column,(n,n))*g))
        end
        orientation = relative_orientation_bases(g_basis1,oriented_basis2)
        # this is the relative orientation of g.first cell compared to second cell
        push!(elements, (g,orientation))
    end
    return elements
end

function boundaries_in_group_ring_with_orientation(cell,basis,cell_dimension,cell_dictionary)
    #= cell dictionary should have the standard cells of the complex, together with a choice of 
        orientation of each such standard cell
    =#
    boundary_cells = []
    lowest_dim = minimum(keys(cell_dictionary))
    #Create a list of the numbers of vertices of cells in codimension 1
    if cell_dimension == lowest_dim
        # only trivial boundaries for cells of lowest dimension
        return boundary_cells
    end
    
    # determine all possible sizes of vertex sets of codim-1 cells
    sizes_in_codimension_1 = Set()
    for (codim1_cell,codim1_basis) in cell_dictionary[cell_dimension-1]
        push!(sizes_in_codimension_1,size(codim1_cell)[2])
    end
    
    vertices = collect(eachcol(cell)) # set of the vertices of cell
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

            # Now find out whether it's in the orbit of a standard cell in codimension 1
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
                        sign = relative_orientation_bases(basis, extended_basis) #this is epsilon(tau',sigma) in Elbaz et al
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


# Stuff that's currently not used:

# pre-check used for speed up by Elbaz et al. 
# Probably won't help much here as there are so few forms and they are almost determined by their determinant.
# Not in use at the moment.
function same_values((form1, min_vectors1),(form2, min_vectors2))
    min_values1 = Set()
    for vector in eachcol(min_vectors1)
        push!(min_values1,transpose(vector)*form1*vector)
    end
    min_values2 = Set()
    for vector in eachcol(min_vectors2)
        push!(min_values2,transpose(vector)*form2*vector)
    end
    return min_values1 == min_values2
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

function a_basis_in_integer_lattices(matrix)
    #= matrix an integer matrix whose columns form a basis of the lattice
    returns one basis
    =#
    columns = collect(eachcol(matrix)) # set of all columns
    dimension = length(columns[1])
    for candidate in Combinatorics.permutations(columns,dimension) # all subsets of columns of with dimension-many elements
        matrix_candidate = transpose(reduce(vcat,transpose.(candidate)))
        if (round(det(matrix_candidate)) == 1) || (round(det(matrix_candidate)) == -1)
            #if that's the case, candidate is a basis of Z^dimension
            return matrix_candidate
        end
    end
end