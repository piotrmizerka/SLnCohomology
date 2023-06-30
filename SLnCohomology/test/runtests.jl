using Groups
using LowCohomologySOS
using PermutationGroups
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
    include("helpful_functions_tests.jl")
    include("permutation_matrices_tests.jl")
end