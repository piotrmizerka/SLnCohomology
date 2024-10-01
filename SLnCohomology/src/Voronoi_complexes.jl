function isposdef_minors(matrix)
    #= isposdef can display incorrect results for singular matrices due to rounding errors. 
        This function is an exact version Sylvester's criterion.
    =#
    @assert (typeof(matrix) == Matrix{Int64} || typeof(matrix) == Matrix{Rational}) 
        "The input needs to be rational for exact computation."
    @assert issymmetric(matrix) "The matrix is not symmetric!"
    for i in 1:size(matrix)[1]
        minor = matrix[1:i, 1:i]
        if !(detx(minor) > 0)
            return false
        end
    end
    return true
end

function add_sl_n_orbits(list_of_forms)
    #= Takes a list of forms, all of the same dimension such that they lie in distinct GL_n orbits. 
    Adds their SL_n orbits ot the list.
    =#
    dim = size(list_of_forms[1])[1]
    flip = flip_matrix(dim)
    if iseven(dim)
        for form in list_of_forms
            SL_orbit = flip*form*flip
            # check whether its in the orbit of one of the other cells
            if isposdef_minors(SL_orbit)
                if ! orbit_in_list(SL_orbit,list_of_forms)
                    push!(list_of_forms,SL_orbit)
                    println("I added an SL_n orbit.")
                end   
            end
        end        
    end
    # if dim is odd, GL and SL orbits are the same, so nothing tbd
    return list_of_forms
end

# Returns the Voronoi cells for SL(n,ℤ)
function cells_sln(n::Integer)
    if n == 2
        A_2 = [2 -1
               -1 2]
        forms = [A_2]
    elseif n == 3
        A_3 = [2 -1 0
               -1 2 -1
               0 -1 2]
        forms = [A_3]
    elseif n == 4
        D_4 = [2 0 1 0
               0 2 -1 0
               1 -1 2 -1
               0 0 -1 2]

        A_4 = [2 -1 0 0
               -1 2 -1 0
               0 -1 2 -1
               0 0 -1 2]

        forms = [D_4, A_4]
    end

    return SLnCohomology.Voronoi_cells(n,forms)
end

function create_polyhedron(cell)
    #= cell is a matrix whose columns are the minimal vectors spanning the corresponding polyhedron
    creates a polyhedron in the polyhedra package
    =#
    # lib = DefaultLibrary{Float64}(GLPK.Optimizer) # one could set different libraries here
    lib = CDDLib.Library(:exact)
    poly_vertices = Vector{Int64}[]
    for col in eachcol(cell)
        push!(poly_vertices,vec(quadratic_form(col)))
    end
    poly_cell = polyhedron(vrep(poly_vertices),lib) 
    return poly_cell
end

function facets(polyhedral_cell, min_vectors, codim_1_cells)
    #=
    min_vectors must be ordered in the same way as the points in polyhedral_cell
    =#
    removehredundancy!(polyhedral_cell) # removes redundant halfspaces s.t. halfspaces become facets
    vertex_list = collect(points(polyhedral_cell))
    for halfspace in eachindex(halfspaces(polyhedral_cell))
        # create a list with all min_vectors that lie on this facet
        vertices_facet = []
        for vertex in incidentpoints(polyhedral_cell, halfspace)
            vertex_index = findfirst(x -> x==vertex, vertex_list)
            push!(vertices_facet,min_vectors[vertex_index])
        end
        # now turn into matrix again to make accessible to other calculations
        facet = transpose(reduce(vcat,transpose.(vertices_facet)))
        if isposdef_minors(quadratic_form(facet))
            # then it intersects the interior non-trivially
            if !orbit_in_list(facet, codim_1_cells)
                push!(codim_1_cells,facet)
            end
        end
    end
    return codim_1_cells
end

function facets_min_vectors_cell(cell,codim_1_cells)
    #= cell is a matrix whose columns are the minimal vectors spanning the corresponding polyhedron
    computes the corresponding polyhedral cell and extracts the minimal vectors in a list
    =#
    return facets(create_polyhedron(cell), collect(eachcol(cell)), codim_1_cells)
end

function flip_matrix(n)
    flip = Matrix(1I, n, n)
    flip[1] = -1
    return flip
end

function Voronoi_cells(n,perfect_forms)
    #= perfect_forms: list of perfect forms in dimension n, up to action of GL
    =#
    dim_symmetric_space = div(n*(n+1),2)-1
    cells_SLn = Dict()
    perfect_forms = add_sl_n_orbits(perfect_forms)
    perfect_forms_min_vec_rep = []
    for perfect_form in perfect_forms
        # compute the minimal vectors and put them in the right format
        push!(perfect_forms_min_vec_rep, transpose(reduce(vcat,transpose.(minimal_vectors(perfect_form)))))
        println("I computed minimal vectors.")
    end
    cells_SLn[dim_symmetric_space] = perfect_forms_min_vec_rep # the perfect forms are the top cells
    for dim in dim_symmetric_space-1:-1:n-1 # these are the dimensions where there are cells that intersect the interior non-trivially
        println("I'm computing in dimension $dim")
        cells_SLn[dim] = []
        for cell in cells_SLn[dim+1]
            facets_min_vectors_cell(cell,cells_SLn[dim])
        end
    end
    return cells_SLn
end