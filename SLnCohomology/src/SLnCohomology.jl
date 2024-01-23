module SLnCohomology

using Combinatorics
using GAP
using Groups
using JSON
using JuMP
using LinearAlgebra
using LowCohomologySOS
using PermutationGroups
using Permutations
using Serialization
using SparseArrays
using Polyhedra

include("helpful_functions.jl")
include("permutation_matrices.jl")
include("Plesken_Souvignier.jl")
include("functions_differential.jl")
include("Voronoi_complexes.jl")

end