N=3 
# parameter N needs to be set globally; if just execute this script, uncomment above

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using Revise
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Serialization
using SLnCohomology

cells_sln = deserialize(joinpath(@__DIR__, "./precomputed_cells/sl"*string(N)*"_cells_new.sjl"))

oriented_cells_sln = SLnCohomology.oriented_cells_dict(cells_sln)

sln_data = Dict()
sln_data["stabilisers"] = SLnCohomology.stabilisers_dict(oriented_cells_sln) # compute the stabilisers
sln_data["boundaries"]  = SLnCohomology.boundaries_dict(oriented_cells_sln) # compute the boundaries

serialize(joinpath(@__DIR__, "precomputed_boundaries/sl"*string(N)*"_bound_stab_new.sjl"), sln_data)