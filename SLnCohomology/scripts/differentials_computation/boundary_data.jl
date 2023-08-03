N=5

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using Revise
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Serialization
using SLnCohomology

cells_sln = deserialize(joinpath(@__DIR__, "./precomputed_cells/sl"*string(N)*"_cells.sjl"))

oriented_cells_sln = Dict()
for dimension in keys(cells_sln)
    oriented_cells_sln[dimension] = []
    for cell in cells_sln[dimension]
        forms_cell = []
        #first convert into vectors (could be combined with below, but easier like this for now)
        for minimal_vector in eachcol(cell)
            push!(forms_cell, vec(SLnCohomology.quadratic_form(minimal_vector)))
        end
        basis = SLnCohomology.extract_basis(forms_cell)
        push!(oriented_cells_sln[dimension],(cell,basis))
    end
end


# Compute all the stabilisers
stabilisers_sln = Dict()
for dimension in keys(cells_sln)
    println("Dimension $dimension")
    cell_list = oriented_cells_sln[dimension]
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
    stabilisers_sln[dimension] = stabilisers_this_dimension
end

# compute the boundaries in SL_N
boundaries_sln = Dict()
for (dimension, cell_list) in oriented_cells_sln
    boundaries_this_dimension = []
    for (cell,basis) in cell_list
        push!(boundaries_this_dimension,SLnCohomology.boundaries_in_group_ring_with_orientation(cell,basis,dimension,oriented_cells_sln))
    end
    boundaries_sln[dimension] = boundaries_this_dimension
end

sln_data = Dict()
sln_data["boundaries"] = boundaries_sln
sln_data["stabilisers"] = stabilisers_sln

serialize(joinpath(@__DIR__, "precomputed_boundaries/sl"*string(N)*"_bound_stab.sjl"), sln_data)