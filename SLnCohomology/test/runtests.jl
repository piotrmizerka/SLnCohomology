using Groups
using LinearAlgebra
using LowCohomologySOS
using Multisets
using PermutationGroups
using Serialization
using SLnCohomology
using SparseArrays
using Test

function cyclic_group(n::Integer)
    A = Alphabet([:a, :A], [2, 1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    return FPGroup(F, [a^n => one(F)])
end

# Save all elts from SL(n,p) and keep their order in the dictionary.
# We use it later in tests to create permutation matrices.
function slnp_ordered(n::Integer,p::Integer)
    slnp_matrices = SLnCohomology.sl_n_p(n,p)
    matrices_orders = Dict(slnp_matrices[i] => i for i in eachindex(slnp_matrices))

    return matrices_orders, slnp_matrices
end

@testset "SLnCohomology" begin
    # include("functions_differential_tests.jl") only include if change to the computation of the chain complex
    include("cohomology_trivial_coefficients_tests.jl")
    include("helpful_functions_tests.jl")
    include("induction_tests.jl")
    include("matrices_to_sln_tests.jl")
    include("subgroup_representations_tests.jl")
    include("Voronoi_complexes_tests.jl")
end