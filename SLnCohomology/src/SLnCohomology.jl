module SLnCohomology

using AbstractAlgebra
using BlockArrays
using Combinatorics
using GAP
using Groups
using JSON
using JuMP
using LinearAlgebra
using LowCohomologySOS
using Permutations
using PermutationGroups
using Polyhedra
using Serialization
using SparseArrays

import GLPK

include("functions_differential.jl")
include("helpful_functions.jl")
include("matrices_to_sln.jl")
include("nontrivial_cohomology.jl")
include("Plesken_Souvignier.jl")
include("Voronoi_complexes.jl")

end