#Pkg.add("Polyhedra")
using Polyhedra

function create_polyhedron(cell)
    poly_vertices = Vector{Int64}[]
    for col in eachcol(cell)
        push!(poly_vertices,vec(SLnCohomology.quadratic_form(col)))
    end
    poly_cell = polyhedron(vrep(poly_vertices))
    return poly_cell
end

function facets(polyhedral_cell, min_vectors, codim_1_cells)
    #=
    min_vectors must be ordered in the same way as the points in polyhedral_cell
    =#
    vertex_list = collect(points(polyhedral_cell))
    for halfspace in eachindex(halfspaces(polyhedral_cell))
        # create a list with all min_vectors that lie on this facet
        vertices_facet = []
        for vertex in incidentpoints(polyhedral_cell, halfspace) # read out using Oscar/polymake in the long term, this might get shorter
            # also careful: At the moment, I'm not reducing the cells, this should be done or avoided by using polymake
            vertex_index = findfirst(x -> x==vertex, vertex_list) # this should somehow get quicker, but I couldn't find out how
            push!(vertices_facet,min_vectors[vertex_index])
        end
        #now turn into matrix again to make accessible to other calculations, orbits extract_basis
        facet = transpose(reduce(vcat,transpose.(vertices_facet)))
        if isposdef(SLnCohomology.quadratic_form(facet))
            # then it intersects the interior non-trivially - make this a separate function to increase readability
            if !orbit_in_list(facet, codim_1_cells)
                push!(codim_1_cells,facet)
            end
        end
    end
    return codim_1_cells
end

function min_vectors_facet(cell,codim_1_cells)
    min_vectors = collect(eachcol(cell))
    return facets(create_polyhedron(cell), min_vectors, codim_1_cells)
end

#=
STILL TBD: The sign flip for even n Implement in the following.
flip_matrix = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

for (dimension,cells) in cells_GL4
    for cell in cells
        SL_orbit = flip_matrix*cell
        # check whether its in the orbit of one of the other cells
        if ! orbit_in_list(SL_orbit,cells_SL4[dimension])
            push!(cells_SL4[dimension],SL_orbit)
        end         
    end
end
=#

function Voronoi_cells(n,perfect_forms)
    #= perfect_forms: list of perfect forms in dimension n, up to action of GL
    =#
    dim_symmetric_space = n*(n+1)/2-1
    cells_SLn = Dict()
    perfect_forms_min_vec_rep = []
    for perfect_form in perfect_forms
        # compute the minimal vectors and put them in the right format
        push!(perfect_forms_min_vec_rep, transpose(reduce(vcat,transpose.(SLnCohomology.minimal_vectors(perfect_form)))))
        println("I computed minimal vectors.")
    end
    cells_SLn[dim_symmetric_space] = perfect_forms_min_vec_rep # the perfect forms are the top cells
    for dim in dim_symmetric_space-1:-1:n-1 # these are the dimensions where there are cells that intersect the interior non-trivially
        println("I'm computing in dimension $dim")
        cells_SLn[dim] = []
        for cell in cells_SLn[dim+1]
            min_vectors_facet(cell,cells_SLn[dim]) # function needs better name
        end
    end
    return cells_SLn
end