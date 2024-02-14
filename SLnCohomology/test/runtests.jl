using Groups
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

@testset "SLnCohomology" begin
    include("functions_differential_tests.jl")
    include("helpful_functions_tests.jl")
    include("matrices_to_sln_tests.jl")
    include("nontrivial_cohomology_tests.jl")
end