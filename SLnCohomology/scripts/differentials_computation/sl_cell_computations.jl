# Compute the cell structure for SL_3, SL_4, SL_5. At the moment, this is actually just the data from Dana Yasaki. Maybe redo on our own at some point.
# It's done very ad hoc at the moment, because it might need to be changed again anyway
# Also important: The 9-dimensional cell for GL_4 is currently missing.


using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using Revise
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Serialization
using SLnCohomology

#SL_3

e1 = [1,0,0]
e2 = [0,1,0]
e3 = [0,0,1]
e12 = [1,1,0]
e13 = [1,0,1]
e123 = [1,1,1]
m2 = [e1 e2 e3] # standard
m3_1 = [e1 e2 e3 e12] # 2-additive
m3_2 = [e1 e2 e3 e123] # 3-additive
m4 = [e1 e2 e3 e12 e13] # double-triple
m5 = [e1 e2 e3 e12 e13 e123] # Q
cells_SL3 = Dict()
cells_SL3[2] = [m2]
cells_SL3[3] = [m3_1,m3_2]
cells_SL3[4] = [m4]
cells_SL3[5] = [m5]

serialize(joinpath(@__DIR__, "precomputed_cells/sl3_cells.sjl"), cells_SL3)

#SL_4

#GL_4 data
# list of vertices
vertices_Voronoi_4 = [[ 0, 0, 0, 1 ],
[ 0, 0, 1, -1 ],
[ 0, 0, 1, 0 ],
[ 0, 1, -1, 0 ],
[ 0, 1, 0, -1 ],
[ 0, 1, 0, 0 ],
[ 0, 1, 1, -1 ],
[ 1, -1, 0, 0 ],
[ 1, 0, -1, 0 ],
[ 1, 0, 0, -1 ],
[ 1, 0, 0, 0 ],
[ 1, 0, 1, -1 ],
[ 1, 1, 0, -1 ]
]

# types of 6-dimensional simplices for GL_4:
index6_1 = [3,5,6,7,8,9,10]
type6_1 = reduce(hcat,[vertices_Voronoi_4[index] for index in index6_1])
index6_2 = [3,5,6,8,9,10,11]
type6_2 = reduce(hcat,[vertices_Voronoi_4[index] for index in index6_2])
index6_3 = [3,5,6,8,10,11,13]
type6_3 = reduce(hcat,[vertices_Voronoi_4[index] for index in index6_3])
index6_4 = [3,5,6,8,9,11,13]
type6_4 = reduce(hcat,[vertices_Voronoi_4[index] for index in index6_4])

# types of 7-dimensional simplices for GL_4:
index7_1 = [3,5,6,7,8,9,10,11]
type7_1 = reduce(hcat,[vertices_Voronoi_4[index] for index in index7_1])
index7_2 = [3,5,6,8,9,10,11,13]
type7_2 = reduce(hcat,[vertices_Voronoi_4[index] for index in index7_2])

# types of 8-dimensional simplices for GL_4:
index8_1 = [3,5,6,7,8,9,10,11,13]
type8_1 = reduce(hcat,[vertices_Voronoi_4[index] for index in index8_1])
index8_2 = [3,6,7,8,9,10,11,12,13]
type8_2 = reduce(hcat,[vertices_Voronoi_4[index] for index in index8_2])

# 9-dimensional simplices missing for now (one of them non-simplicial)

cells_GL4 = Dict()
cells_GL4[6] = [type6_1,type6_2,type6_3,type6_4]
cells_GL4[7] = [type7_1,type7_2]
cells_GL4[8] = [type8_1,type8_2]

cells_SL4 = Dict()
cells_SL4[6] = [type6_1,type6_2,type6_3,type6_4]
cells_SL4[7] = [type7_1,type7_2]
cells_SL4[8] = [type8_1,type8_2]

flip_matrix = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
# multiplying with this matrix gives the potential second SL_n orbit contained in a GL_n orbit
# though this is somewhat redundant here because "by chance" the SL_n orbits seem to be the same
# keep it in the code in case this needs to be redone for higher n

function same_orbit(matrix1,matrix2)
    #= Checks whether matrix1 lies in the SL_n orbit as matrix2
    =#
    if length(SLnCohomology.stabiliser_coset_SL(SLnCohomology.quadratic_form(matrix1),SLnCohomology.quadratic_form(matrix2))) > 0
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


for (dimension,cells) in cells_GL4
    for cell in cells
        SL_orbit = flip_matrix*cell
        # check whether its in the orbit of one of the other cells
        if ! orbit_in_list(SL_orbit,cells_SL4[dimension])
            push!(cells_SL4[dimension],SL_orbit)
        end         
    end
end

#compute boundaries of these SL cells
# only works for SL_4 because everything is simplicial, would need to be redone for higher even n

cell_orbits_in_dimension_5 = Matrix{Int64}[]
for higher_cell in cells_SL4[6]
    for j in 1:size(higher_cell)[2]
        boundary = higher_cell[begin:end,begin:end .!=j] # exclude the j-th column
        if rank(boundary) == 4
            # only go on if it lies in the interior of the symmetric space
            if !orbit_in_list(boundary,cell_orbits_in_dimension_5)
                push!(cell_orbits_in_dimension_5,boundary)
            end
        end
    end
end
cells_SL4[5] = cell_orbits_in_dimension_5
dimension = 5
number_cells = length(cell_orbits_in_dimension_5)
println("Number of cells in dimension $dimension: $number_cells")

cell_orbits_in_dimension_4 = Matrix{Int64}[]
for higher_cell in cells_SL4[5]
    for j in 1:size(higher_cell)[2]
        boundary = higher_cell[begin:end,begin:end .!=j] # exclude the j-th column
        if rank(boundary) == 4
            # only go on if it lies in the interior of the symmetric space
            if !orbit_in_list(boundary,cell_orbits_in_dimension_4)
                push!(cell_orbits_in_dimension_4,boundary)
            end
        end
    end
end
cells_SL4[4] = cell_orbits_in_dimension_4
dimension = 4
number_cells = length(cell_orbits_in_dimension_4)
println("Number of cells in dimension $dimension: $number_cells")

cell_orbits_in_dimension_3 = Matrix{Int64}[]
for higher_cell in cells_SL4[4]
    for j in 1:size(higher_cell)[2]
        boundary = higher_cell[begin:end,begin:end .!=j] # exclude the j-th column
        if rank(boundary) == 4
            # only go on if it lies in the interior of the symmetric space
            if !orbit_in_list(boundary,cell_orbits_in_dimension_3)
                push!(cell_orbits_in_dimension_3,boundary)
            end
        end
    end
end
cells_SL4[3] = cell_orbits_in_dimension_3

dimension = 3
number_cells = length(cell_orbits_in_dimension_3)
println("Number of cells in dimension $dimension: $number_cells")

serialize(joinpath(@__DIR__, "precomputed_cells/sl4_cells.sjl"), cells_SL4)

#SL_5

# SL_5 data
# list of vertices
vertices_Voronoi_5 = [
[ 0, 0, 0, 0, 1 ],
[ 0, 0, 0, 1, -1 ],
[ 0, 0, 0, 1, 0 ],
[ 0, 0, 1, -1, 0 ],
[ 0, 0, 1, 0, -1 ],
[ 0, 0, 1, 0, 0 ],
[ 0, 0, 1, 0, 1 ],
[ 0, 0, 1, 1, -1 ],
[ 0, 0, 1, 1, 0 ],
[ 0, 1, -1, 0, 0 ],
[ 0, 1, 0, -1, 0 ],
[ 0, 1, 0, 0, -1 ],
[ 0, 1, 0, 0, 0 ],
[ 0, 1, 0, 0, 1 ],
[ 0, 1, 0, 1, -1 ],
[ 0, 1, 0, 1, 0 ],
[ 0, 1, 1, 0, -1 ],
[ 0, 1, 1, 1, 1 ],
[ 1, -1, 0, 0, 0 ],
[ 1, 0, -1, 0, 0 ],
[ 1, 0, 0, -1, 0 ],
[ 1, 0, 0, 0, -1 ],
[ 1, 0, 0, 0, 0 ],
[ 1, 0, 0, 0, 1 ],
[ 1, 0, 0, 1, -1 ],
[ 1, 0, 0, 1, 0 ],
[ 1, 0, 1, 0, -1 ],
[ 1, 0, 1, 1, 1 ],
[ 1, 1, 0, 0, -1 ],
[ 1, 1, 0, 1, 1 ],
[ 1, 1, 1, 1, 1 ]
]

# the following is a list containing a representative of each orbit of cells
# in the tesselation. It starts with cells of dimension 4 and is ordered
# by dimension.
# ATTENTION: The indexing here starts at 1, not at 0.
indices_k_cells = [
[
[ 3, 6, 8, 12, 19 ],
[ 8, 12, 13, 20, 21 ]
],
[
[ 3, 6, 8, 12, 13, 19 ],
[ 3, 8, 12, 13, 15, 19 ],
[ 3, 6, 8, 13, 15, 19 ],
[ 6, 8, 12, 13, 20, 21 ],
[ 3, 8, 13, 15, 20, 21 ]
],
[
[ 3, 6, 8, 12, 13, 15, 19 ],
[ 3, 6, 8, 12, 13, 19, 20 ],
[ 3, 6, 8, 13, 15, 19, 20 ],
[ 3, 6, 12, 13, 15, 19, 20 ],
[ 3, 6, 8, 12, 13, 20, 21 ],
[ 3, 6, 8, 13, 15, 20, 21 ],
[ 3, 8, 12, 13, 15, 20, 21 ],
[ 3, 6, 8, 12, 15, 19, 21 ],
[ 3, 6, 8, 12, 15, 19, 22 ],
[ 3, 6, 12, 19, 21, 22, 23 ]
],
[
[ 3, 6, 8, 12, 13, 15, 19, 20 ],
[ 3, 6, 8, 12, 13, 19, 20, 21 ],
[ 3, 6, 8, 12, 13, 15, 20, 21 ],
[ 3, 6, 8, 12, 13, 15, 19, 21 ],
[ 3, 6, 8, 13, 15, 19, 20, 21 ],
[ 3, 6, 8, 12, 13, 15, 19, 22 ],
[ 3, 6, 12, 13, 15, 19, 21, 22 ],
[ 3, 6, 8, 12, 19, 21, 22, 23 ],
[ 3, 6, 8, 12, 19, 20, 21, 23 ],
[ 3, 6, 12, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 20, 21, 23 ],
[ 3, 6, 12, 13, 20, 21, 22, 23 ],
[ 3, 6, 12, 13, 19, 21, 22, 23 ],
[ 3, 6, 8, 13, 19, 20, 21, 23 ],
[ 3, 6, 12, 13, 19, 20, 21, 29 ],
[ 3, 6, 12, 13, 19, 22, 23, 29 ]
],
[
[ 3, 6, 8, 12, 13, 15, 19, 20, 21 ],
[ 3, 6, 8, 12, 13, 15, 19, 21, 22 ],
[ 3, 8, 12, 13, 15, 19, 20, 21, 22 ],
[ 3, 6, 8, 13, 15, 19, 20, 21, 22 ],
[ 3, 6, 8, 12, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 19, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 19, 20, 21, 23 ],
[ 3, 6, 12, 13, 19, 20, 21, 22, 23 ],
[ 6, 8, 12, 13, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 13, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 15, 20, 22, 23 ],
[ 3, 6, 12, 13, 15, 20, 21, 22, 23 ],
[ 3, 6, 12, 13, 15, 19, 21, 22, 23 ],
[ 3, 6, 12, 13, 15, 19, 20, 22, 23 ],
[ 3, 6, 8, 12, 13, 19, 20, 21, 29 ],
[ 3, 6, 12, 13, 19, 20, 21, 22, 29 ],
[ 3, 6, 12, 13, 19, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 13, 19, 22, 23, 29 ],
[ 6, 8, 12, 13, 19, 21, 22, 23, 29 ],
[ 3, 6, 12, 13, 19, 20, 21, 23, 29 ],
[ 3, 6, 12, 13, 15, 20, 22, 23, 29 ],
[ 2, 6, 8, 12, 13, 19, 22, 23, 29 ]
],
[
[ 3, 6, 8, 12, 13, 15, 19, 20, 21, 22 ],
[ 3, 6, 8, 12, 13, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 15, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 15, 19, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 15, 19, 20, 22, 23 ],
[ 3, 6, 12, 13, 15, 19, 20, 21, 22, 23 ],
[ 3, 8, 12, 13, 15, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 13, 15, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 19, 20, 21, 22, 29 ],
[ 3, 6, 8, 12, 13, 19, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 13, 19, 20, 21, 23, 29 ],
[ 3, 6, 12, 13, 19, 20, 21, 22, 23, 29 ],
[ 3, 6, 12, 13, 15, 19, 20, 21, 22, 29 ],
[ 3, 8, 12, 13, 15, 19, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 13, 15, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 13, 15, 20, 22, 23, 29 ],
[ 3, 6, 12, 13, 15, 19, 20, 22, 23, 29 ],
[ 3, 8, 13, 15, 19, 21, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 13, 19, 21, 22, 23, 29 ],
[ 2, 6, 8, 12, 13, 15, 19, 21, 23, 29 ],
[ 2, 6, 8, 12, 13, 15, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 13, 19, 20, 22, 23, 25 ],
[ 3, 6, 8, 12, 19, 20, 22, 23, 25, 29 ],
[ 6, 8, 13, 15, 20, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 13, 15, 17, 19, 20, 21, 22 ]
],
[
[ 3, 6, 8, 12, 13, 15, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 13, 19, 20, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 13, 15, 19, 20, 21, 22, 29 ],
[ 3, 6, 8, 12, 13, 15, 19, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 13, 15, 19, 20, 22, 23, 29 ],
[ 3, 6, 12, 13, 15, 19, 20, 21, 22, 23, 29 ],
[ 3, 6, 8, 13, 15, 19, 21, 22, 23, 25, 29 ],
[ 2, 3, 6, 11, 12, 13, 15, 19, 21, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 13, 15, 19, 21, 22, 23, 29 ],
[ 2, 6, 8, 12, 15, 19, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 12, 13, 15, 19, 20, 22, 23, 25 ],
[ 3, 6, 8, 12, 13, 19, 20, 22, 23, 25, 29 ],
[ 3, 6, 8, 12, 13, 15, 19, 20, 23, 25, 29 ],
[ 3, 6, 8, 12, 13, 15, 20, 22, 23, 25, 29 ],
[ 3, 6, 8, 13, 15, 19, 20, 21, 22, 23, 25 ],
[ 3, 6, 8, 13, 15, 20, 21, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 13, 19, 20, 22, 23, 25, 29 ],
[ 3, 6, 8, 12, 13, 15, 17, 19, 20, 21, 22 ],
[ 3, 6, 8, 12, 15, 17, 19, 20, 21, 22, 25 ],
[ 6, 8, 12, 15, 17, 19, 20, 21, 22, 23, 25 ],
[ 3, 6, 8, 15, 17, 19, 20, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 17, 19, 20, 21, 23, 25, 29 ],
[ 2, 6, 8, 11, 12, 17, 19, 20, 23, 25, 29 ]
],
[
[ 3, 6, 8, 12, 13, 15, 19, 20, 21, 22, 23, 29 ],
[ 2, 3, 6, 8, 11, 12, 13, 15, 19, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 12, 13, 15, 19, 20, 22, 23, 25, 29 ],
[ 3, 6, 8, 13, 15, 19, 20, 21, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 13, 15, 19, 20, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 15, 19, 20, 21, 22, 23, 25, 29 ],
[ 2, 3, 6, 8, 11, 12, 20, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 11, 12, 13, 15, 19, 20, 21, 22, 29 ],
[ 3, 6, 8, 12, 13, 15, 17, 19, 20, 21, 22, 23 ],
[ 3, 6, 8, 12, 15, 17, 19, 20, 21, 22, 23, 25 ],
[ 3, 6, 8, 15, 17, 19, 20, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 12, 15, 17, 19, 20, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 15, 17, 19, 20, 21, 22, 25, 29 ],
[ 2, 6, 8, 12, 15, 17, 19, 20, 21, 23, 25, 29 ],
[ 2, 6, 8, 11, 12, 17, 19, 20, 21, 23, 25, 29 ],
[ 3, 6, 15, 17, 19, 20, 21, 22, 23, 25, 27, 29 ]
],
[
[ 2, 3, 6, 8, 11, 12, 13, 15, 19, 20, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 12, 13, 15, 17, 19, 20, 21, 22, 23, 29 ],
[ 3, 6, 8, 12, 15, 17, 19, 20, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 13, 15, 17, 19, 20, 21, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 15, 17, 19, 20, 21, 22, 23, 25, 29 ],
[ 2, 6, 8, 12, 13, 15, 17, 19, 20, 22, 23, 25, 29 ],
[ 2, 6, 8, 11, 12, 17, 19, 20, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 15, 17, 19, 20, 21, 22, 23, 25, 27, 29 ],
[ 3, 6, 12, 15, 17, 19, 20, 21, 22, 23, 25, 27, 29 ]
],
[
[ 2, 3, 6, 8, 11, 12, 13, 15, 17, 19, 20, 21, 22, 23, 25, 29 ],
[ 3, 6, 8, 12, 15, 17, 19, 20, 21, 22, 23, 25, 27, 29 ],
[ 3, 6, 8, 13, 15, 17, 19, 20, 21, 22, 23, 25, 27, 29 ],
[ 3, 4, 6, 12, 15, 17, 19, 20, 21, 22, 23, 25, 27, 29 ]
],
[
[ 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 15, 17, 19, 20, 21, 22, 23, 25, 27, 29 ],
[ 1, 3, 6, 7, 9, 13, 14, 16, 18, 23, 24, 26, 28, 30, 31 ],
[ 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 19, 20, 21, 22, 23 ]
]
]

cells_GL5 = Dict()
for dimension in range(4,14)
    cells_current_dimension = []
    for cell_indices in indices_k_cells[dimension-3]
        cell = reduce(hcat,[vertices_Voronoi_5[index] for index in cell_indices])
        push!(cells_current_dimension, cell)
    end
    cells_GL5[dimension] = cells_current_dimension
    number_cells = length(cells_current_dimension)
    println("Number of cells in dimension $dimension: $number_cells")
end

#odd dimension: SL orbits are the same
cells_SL5 = cells_GL5

serialize(joinpath(@__DIR__, "precomputed_cells/sl5_cells.sjl"), cells_SL5)