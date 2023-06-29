include("../../src/permutation_matrices.jl")

using LinearAlgebra
using IterTools
using PermutationGroups
using SparseArrays

const N = 3
const p = 3

slnp_dict, matrices_fixed_order = slnp(N,p)

function permutation_matrices(matrices)
    result = Dict()
    it = 0
    for M in matrices
        left_action_matrix(i) = projection(M*matrices[i],p)
        perm = PermutationGroups.Perm([slnp_dict[left_action_matrix(i)] for i in eachindex(matrices)])
        result[M] = permutation_matrix(perm)
        if it%1000 == 0
            @info it
        end
        it += 1
    end
    return result
end

perm_mats = permutation_matrices(matrices_fixed_order)
deg = slnp_order(N,p)

# Permutation matrix corresponding do a group ring elt ξ given by the perm repr
function representing_matrix(ξ,p::Integer)
    coeff_type = typeof(first(ξ.coeffs))
    result = sparse(coeff_type.(zeros(deg,deg)))
    RG = parent(ξ)
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ.coeffs[i]*perm_mats[proj]
    end
    return result
end

# For SL(3,Z), modular projection with p = 3 ##############################################################################
include("../differentials_computation/sl3_laplacians.jl");
@info "SL(3,Z) via SL(3,3):"
# A sanity check - since H^1(SL(3,Z),π) = 0, we shall get the full rank for the perm repr of Δ₄:
delta_4_perm = Array(representing_matrix(Δ₄[1,1],p))
@assert delta_4_perm' == delta_4_perm
@info "rank(H^1) = ", size(delta_4_perm)[1]-rank(delta_4_perm)

# A sanity check: since the permutation representation has a nontrivial fixed pont set, 
# we shall get nontrivial zero cohomology for this representation.
delta_5_perm = Array(representing_matrix(Δ₅[1,1],p))
@assert delta_5_perm' == delta_5_perm
@info "rank(H^0) = ", size(delta_5_perm)[1]-rank(delta_5_perm) # non-full rank - means nontrivial 0-cohomology

# H^2 - we get nontrivial cohomology!
delta_3_11 = Array(representing_matrix(Δ₃[1,1],p))
delta_3_12 = Array(representing_matrix(Δ₃[1,2],p))
delta_3_21 = Array(representing_matrix(Δ₃[2,1],p))
delta_3_22 = Array(representing_matrix(Δ₃[2,2],p))
delta_3_perm = [hcat(delta_3_11,delta_3_12);hcat(delta_3_21,delta_3_22)]
@assert delta_3_perm' == delta_3_perm
@info "rank(H^2) = ", size(delta_3_perm)[1]-rank(delta_3_perm) # non-full rank - means nontrivial 2-cohomology
##########################################################################################
