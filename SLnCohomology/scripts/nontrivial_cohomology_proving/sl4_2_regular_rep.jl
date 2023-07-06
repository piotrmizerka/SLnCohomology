include("../../src/permutation_matrices.jl")

using LinearAlgebra
using PermutationGroups
using SLnCohomology
using SparseArrays

const N = 4
const p = 2

slnp_dict, matrices_fixed_order = slnp(N,p)
perm_mats = SLnCohomology.permutation_matrices(matrices_fixed_order, slnp_dict, p)

# For SL(4,Z), modular projection with p = 2 ##############################################################################
function is_singular(U) # U is an upper-triangular square matrix
    entry_type = typeof(first(U))
    for i in 1:size(U)[1]
        if U[i,i] == zero(entry_type)
            return true
        end
    end
    return false
end

include("../differentials_computation/sl4_laplacians.jl");
@info "SL(4,Z) via SL(4,2):"
# A sanity check - since, by Bader-Sauer, H^2(SL(4,Z),π) = 0,
#  we shall get the full rank for the perm repr of Δ₄:
delta_7_11 = SLnCohomology.representing_matrix(Δ₇[1,1],p,perm_mats)
delta_7_12 = SLnCohomology.representing_matrix(Δ₇[1,2],p,perm_mats)
delta_7_21 = SLnCohomology.representing_matrix(Δ₇[2,1],p,perm_mats)
delta_7_22 = SLnCohomology.representing_matrix(Δ₇[2,2],p,perm_mats)
delta_7_perm = [hcat(delta_7_11,delta_7_12);hcat(delta_7_21,delta_7_22)]
@assert delta_7_perm' == delta_7_perm
@info "Is H^2 nontrivial?", is_singular(qr(Float64.(delta_7_perm)).R)

# H^3 - we get nontrivial cohomology!
delta_6_11 = SLnCohomology.representing_matrix(Δ₆[1,1],p,perm_mats)
delta_6_12 = SLnCohomology.representing_matrix(Δ₆[1,2],p,perm_mats)
delta_6_13 = SLnCohomology.representing_matrix(Δ₆[1,3],p,perm_mats)
delta_6_14 = SLnCohomology.representing_matrix(Δ₆[1,4],p,perm_mats)
delta_6_21 = SLnCohomology.representing_matrix(Δ₆[2,1],p,perm_mats)
delta_6_22 = SLnCohomology.representing_matrix(Δ₆[2,2],p,perm_mats)
delta_6_23 = SLnCohomology.representing_matrix(Δ₆[2,3],p,perm_mats)
delta_6_24 = SLnCohomology.representing_matrix(Δ₆[2,4],p,perm_mats)
delta_6_31 = SLnCohomology.representing_matrix(Δ₆[3,1],p,perm_mats)
delta_6_32 = SLnCohomology.representing_matrix(Δ₆[3,2],p,perm_mats)
delta_6_33 = SLnCohomology.representing_matrix(Δ₆[3,3],p,perm_mats)
delta_6_34 = SLnCohomology.representing_matrix(Δ₆[3,4],p,perm_mats)
delta_6_41 = SLnCohomology.representing_matrix(Δ₆[4,1],p,perm_mats)
delta_6_42 = SLnCohomology.representing_matrix(Δ₆[4,2],p,perm_mats)
delta_6_43 = SLnCohomology.representing_matrix(Δ₆[4,3],p,perm_mats)
delta_6_44 = SLnCohomology.representing_matrix(Δ₆[4,4],p,perm_mats)
delta_6_perm = [
    hcat(delta_6_11,delta_6_12,delta_6_13,delta_6_14);
    hcat(delta_6_21,delta_6_22,delta_6_23,delta_6_24);
    hcat(delta_6_31,delta_6_32,delta_6_33,delta_6_34);
    hcat(delta_6_41,delta_6_42,delta_6_43,delta_6_44)
]
@assert delta_6_perm' == delta_6_perm
@info "Is H^3 nontrivial?", is_singular(qr(Float64.(delta_6_perm)).R)

# H^4 - here our perm rep yields trivial cohomology
delta_5_11 = SLnCohomology.representing_matrix(Δ₅[1,1],p,perm_mats)
delta_5_12 = SLnCohomology.representing_matrix(Δ₅[1,2],p,perm_mats)
delta_5_13 = SLnCohomology.representing_matrix(Δ₅[1,3],p,perm_mats)
delta_5_14 = SLnCohomology.representing_matrix(Δ₅[1,4],p,perm_mats)
delta_5_21 = SLnCohomology.representing_matrix(Δ₅[2,1],p,perm_mats)
delta_5_22 = SLnCohomology.representing_matrix(Δ₅[2,2],p,perm_mats)
delta_5_23 = SLnCohomology.representing_matrix(Δ₅[2,3],p,perm_mats)
delta_5_24 = SLnCohomology.representing_matrix(Δ₅[2,4],p,perm_mats)
delta_5_31 = SLnCohomology.representing_matrix(Δ₅[3,1],p,perm_mats)
delta_5_32 = SLnCohomology.representing_matrix(Δ₅[3,2],p,perm_mats)
delta_5_33 = SLnCohomology.representing_matrix(Δ₅[3,3],p,perm_mats)
delta_5_34 = SLnCohomology.representing_matrix(Δ₅[3,4],p,perm_mats)
delta_5_41 = SLnCohomology.representing_matrix(Δ₅[4,1],p,perm_mats)
delta_5_42 = SLnCohomology.representing_matrix(Δ₅[4,2],p,perm_mats)
delta_5_43 = SLnCohomology.representing_matrix(Δ₅[4,3],p,perm_mats)
delta_5_44 = SLnCohomology.representing_matrix(Δ₅[4,4],p,perm_mats)
delta_5_perm = [
    hcat(delta_5_11,delta_5_12,delta_5_13,delta_5_14);
    hcat(delta_5_21,delta_5_22,delta_5_23,delta_5_24);
    hcat(delta_5_31,delta_5_32,delta_5_33,delta_5_34);
    hcat(delta_5_41,delta_5_42,delta_5_43,delta_5_44)
]
# @assert delta_5_perm' == delta_5_perm
@info "Is H^4 nontrivial?", is_singular(qr(Float64.(delta_5_perm)).R)

# H^5:
delta_4_11 = SLnCohomology.representing_matrix(Δ₄[1,1],p,perm_mats)
delta_4_12 = SLnCohomology.representing_matrix(Δ₄[1,2],p,perm_mats)
delta_4_13 = SLnCohomology.representing_matrix(Δ₄[1,3],p,perm_mats)
delta_4_21 = SLnCohomology.representing_matrix(Δ₄[2,1],p,perm_mats)
delta_4_22 = SLnCohomology.representing_matrix(Δ₄[2,2],p,perm_mats)
delta_4_23 = SLnCohomology.representing_matrix(Δ₄[2,3],p,perm_mats)
delta_4_31 = SLnCohomology.representing_matrix(Δ₄[3,1],p,perm_mats)
delta_4_32 = SLnCohomology.representing_matrix(Δ₄[3,2],p,perm_mats)
delta_4_33 = SLnCohomology.representing_matrix(Δ₄[3,3],p,perm_mats)
delta_4_perm = [
    hcat(delta_4_11,delta_4_12,delta_4_13);
    hcat(delta_4_21,delta_4_22,delta_4_23);
    hcat(delta_4_31,delta_4_32,delta_4_33)
]
# @assert delta_4_perm' == delta_4_perm
@info "Is H^5 nontrivial?", is_singular(qr(Float64.(delta_4_perm)).R)
##########################################################################################
