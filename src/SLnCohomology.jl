module SLnCohomology

using AbstractAlgebra
using BlockArrays
using CDDLib
using Combinatorics
using GAP
using Groups
using IntervalLinearAlgebra
using JSON
using JuMP
using LinearAlgebra
using LinearAlgebraX
using LowCohomologySOS
using Permutations
using PermutationGroups
using Polyhedra
using Serialization
using SparseArrays

import GLPK

include("functions_differential.jl")
include("helpful_functions.jl")
include("induction.jl")
include("matrices_to_sln.jl")
include("Plesken_Souvignier.jl")
include("subgroup_representations.jl")
include("Voronoi_complexes.jl")

end