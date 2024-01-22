using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using Revise
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Serialization
using SLnCohomology

#SL_3
A_3 = [2 -1 0
-1 2 -1
0 -1 2]
forms_3 = [A_3]

#SL_4
D_4 = [2 0 1 0
0 2 -1 0
1 -1 2 -1
0 0 -1 2]

A_4 = [2 -1 0 0
-1 2 -1 0
0 -1 2 -1
0 0 -1 2]

forms_4 = [D_4, A_4]

#SL_5
D_5 = [2 0 1 0 0
0 2 -1 0 0
1 -1 2 -1 0
0 0 -1 2 -1
0 0 0 -1 2]

A_5_plus3 = [6 -3 0 0 0
-3 6 -3 0 3
0 -3 6 -3 0
0 0 -3 6 0
0 3 0 0 4]

A_5 = [2 -1 0 0 0
-1 2 -1 0 0
0 -1 2 -1 0
0 0 -1 2 -1
0 0 0 -1 2]

forms_5 = [D_5, A_5_plus3, A_5]

cells_SL3 = SLnCohomology.Voronoi_cells(3,forms_3)
serialize(joinpath(@__DIR__, "precomputed_cells/sl3_cells_new.sjl"), cells_SL3)

cells_SL4 = SLnCohomology.Voronoi_cells(4,forms_4)
serialize(joinpath(@__DIR__, "precomputed_cells/sl4_cells_new.sjl"), cells_SL4)

cells_SL5 = SLnCohomology.Voronoi_cells(5,forms_5)
serialize(joinpath(@__DIR__, "precomputed_cells/sl5_cells_new.sjl"), cells_SL5)