include("../../src/permutation_matrices.jl")

using LinearAlgebra
using IterTools
using PermutationGroups
using SLnCohomology
using SparseArrays

const N = 3
const p = 3

slnp_dict, matrices_fixed_order = slnp(N,p)
perm_mats = SLnCohomology.permutation_matrices(matrices_fixed_order, slnp_dict, p)

# For SL(3,Z), modular projection with p = 3 ##############################################################################
include("../differentials_computation/sl3_laplacians.jl");
@info "SL(3,Z) via SL(3,3):"
# A sanity check - since H^1(SL(3,Z),π) = 0, we shall get the full rank for the perm repr of Δ₄:
delta_4_perm = Array(SLnCohomology.representing_matrix(Δ₄[1,1],p,perm_mats))
@assert delta_4_perm' == delta_4_perm
@info "rank(H^1) = ", size(delta_4_perm)[1]-rank(delta_4_perm)

# A sanity check: since the permutation representation has a nontrivial fixed pont set, 
# we shall get nontrivial zero cohomology for this representation.
delta_5_perm = Array(SLnCohomology.representing_matrix(Δ₅[1,1],p,perm_mats))
@assert delta_5_perm' == delta_5_perm
@info "rank(H^0) = ", size(delta_5_perm)[1]-rank(delta_5_perm) # non-full rank - means nontrivial 0-cohomology

# H^2 - we get nontrivial cohomology!
delta_3_11 = Array(SLnCohomology.representing_matrix(Δ₃[1,1],p,perm_mats))
delta_3_12 = Array(SLnCohomology.representing_matrix(Δ₃[1,2],p,perm_mats))
delta_3_21 = Array(SLnCohomology.representing_matrix(Δ₃[2,1],p,perm_mats))
delta_3_22 = Array(SLnCohomology.representing_matrix(Δ₃[2,2],p,perm_mats))
delta_3_perm = [hcat(delta_3_11,delta_3_12);hcat(delta_3_21,delta_3_22)]
@assert delta_3_perm' == delta_3_perm
@info "rank(H^2) = ", size(delta_3_perm)[1]-rank(delta_3_perm) # non-full rank - means nontrivial 2-cohomology
##########################################################################################
