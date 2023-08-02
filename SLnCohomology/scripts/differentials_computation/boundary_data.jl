N=5

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using Revise
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Serialization
using SLnCohomology

deserialize(joinpath(@__DIR__, "./precomputed_cells/sl"*string(N)*"_cells.sjl"))

# as 5 is odd, we don't need to distinguish between SL5 and Gl5 orbits
# maybe dependent on n later on
cells_SL5 = cells_GL5

oriented_cells_SL5 = Dict()
for dimension in range(4,14)
    oriented_cells_SL5[dimension] = []
    for cell in cells_SL5[dimension]
        forms_cell = []
        #first convert into vectors (could be combined with below, but easier like this for now)
        for minimal_vector in eachcol(cell)
            push!(forms_cell, vec(SLnCohomology.quadratic_form(minimal_vector)))
        end
        basis = SLnCohomology.extract_basis(forms_cell)
        push!(oriented_cells_SL5[dimension],(cell,basis))
    end
end


# Compute all the stabilisers
stabilisers_SL5 = Dict()
for dimension in range(4,14)
    println("Dimension $dimension")
    cell_list = oriented_cells_SL5[dimension]
    stabilisers_this_dimension = []
    cell_count = 0
    for (cell,basis) in cell_list
        cell_count += 1
        println("Cell number $cell_count")
        cell_stabiliser = []
        for (g,orientation) in SLnCohomology.stabiliser_coset_with_orientation((cell,basis), (cell,basis))
            push!(cell_stabiliser,(g,orientation))
        end
        push!(stabilisers_this_dimension,cell_stabiliser)
    end
    stabilisers_SL5[dimension] = stabilisers_this_dimension
end

# compute the boundaries in SL5
boundaries_SL5 = Dict()
for (dimension, cell_list) in oriented_cells_SL5
    boundaries_this_dimension = []
    for (cell,basis) in cell_list
        push!(boundaries_this_dimension,SLnCohomology.boundaries_in_group_ring_with_orientation(cell,basis,dimension,oriented_cells_SL5))
    end
    boundaries_SL5[dimension] = boundaries_this_dimension
end

sl5_data = Dict()
sl5_data["boundaries"] = boundaries_SL5
sl5_data["stabilisers"] = stabilisers_SL5

serialize(joinpath(@__DIR__, "precomputed_boundaries/sl5_bound_stab.sjl"), sl5_data)