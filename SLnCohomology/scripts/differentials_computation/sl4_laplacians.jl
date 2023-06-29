# Let's load the precompuded differential data (can take about 10mins)
include("sl4_utils.jl");

# Sanity checks for vanishing compositions of differentials
d₄x = LowCohomologySOS.embed.(identity, d₄, Ref(RG_Δ[4]));
d₅x = LowCohomologySOS.embed.(identity, d₅, Ref(RG_Δ[4]));
@info d₄x*d₅x == [zero(RG_Δ[4]) zero(RG_Δ[4]) zero(RG_Δ[4]) zero(RG_Δ[4])]
d₅x = LowCohomologySOS.embed.(identity, d₅, Ref(RG_Δ[5]));
d₆x = LowCohomologySOS.embed.(identity, d₆, Ref(RG_Δ[5]));
@info d₅x*d₆x == [
    zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]);
    zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]);
    zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5])
]
d₆x = LowCohomologySOS.embed.(identity, d₆, Ref(RG_Δ[6]));
d₇x = LowCohomologySOS.embed.(identity, d₇, Ref(RG_Δ[6]));
@info d₆x*d₇x == [
    zero(RG_Δ[6]) zero(RG_Δ[6]);
    zero(RG_Δ[6]) zero(RG_Δ[6]);
    zero(RG_Δ[6]) zero(RG_Δ[6]);
    zero(RG_Δ[6]) zero(RG_Δ[6])
]
d₇x = LowCohomologySOS.embed.(identity, d₇, Ref(RG_Δ[7]));
d₈x = LowCohomologySOS.embed.(identity, d₈, Ref(RG_Δ[7]));
@info d₇x*d₈x == [
    zero(RG_Δ[7]) zero(RG_Δ[7]);
    zero(RG_Δ[7]) zero(RG_Δ[7]);
    zero(RG_Δ[7]) zero(RG_Δ[7]);
    zero(RG_Δ[7]) zero(RG_Δ[7])
]

# The stabiliser parts which we have to add to get free modules.
# Hopefully the stabilisers' elements belong to the half_bases.
dim3_stab_part = reshape([one(RG_Δ[3])-averaged_rep(m3_stabs[1], half_basis_Δ[3], RG_Δ[3])], 1, 1)
dim4_stab_part = [
    one(RG_Δ[4])-averaged_rep(m4_stabs[1], half_basis_Δ[4], RG_Δ[4]) zero(RG_Δ[4]) zero(RG_Δ[4]);
    zero(RG_Δ[4]) one(RG_Δ[4])-averaged_rep(m4_stabs[2], half_basis_Δ[4], RG_Δ[4]) zero(RG_Δ[4]);
    zero(RG_Δ[4]) zero(RG_Δ[4]) one(RG_Δ[4])-averaged_rep(m4_stabs[3], half_basis_Δ[4], RG_Δ[4]);
]
dim5_stab_part = [
    one(RG_Δ[5])-averaged_rep(m5_stabs[1], half_basis_Δ[5], RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]);
    zero(RG_Δ[5]) one(RG_Δ[5])-averaged_rep(m5_stabs[2], half_basis_Δ[5], RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]);
    zero(RG_Δ[5]) zero(RG_Δ[5]) one(RG_Δ[5])-averaged_rep(m5_stabs[3], half_basis_Δ[5], RG_Δ[5]) zero(RG_Δ[5]);
    zero(RG_Δ[5]) zero(RG_Δ[5]) zero(RG_Δ[5]) one(RG_Δ[5])-averaged_rep(m5_stabs[4], half_basis_Δ[5], RG_Δ[5])
]
dim6_stab_part = [
    one(RG_Δ[6])-averaged_rep(m6_stabs[1], half_basis_Δ[6], RG_Δ[6]) zero(RG_Δ[6]) zero(RG_Δ[6]) zero(RG_Δ[6]);
    zero(RG_Δ[6]) one(RG_Δ[6])-averaged_rep(m6_stabs[2], half_basis_Δ[6], RG_Δ[6]) zero(RG_Δ[6]) zero(RG_Δ[6]);
    zero(RG_Δ[6]) zero(RG_Δ[6]) one(RG_Δ[6])-averaged_rep(m6_stabs[3], half_basis_Δ[6], RG_Δ[6]) zero(RG_Δ[6]);
    zero(RG_Δ[6]) zero(RG_Δ[6]) zero(RG_Δ[6]) one(RG_Δ[6])-averaged_rep(m6_stabs[4], half_basis_Δ[6], RG_Δ[6])
]
dim7_stab_part = [
    one(RG_Δ[7])-averaged_rep(m7_stabs[1], half_basis_Δ[7], RG_Δ[7]) zero(RG_Δ[7]);
    zero(RG_Δ[7]) one(RG_Δ[7])-averaged_rep(m7_stabs[2], half_basis_Δ[7], RG_Δ[7])
];

# Compute the Laplacians (can take about 30mins).
# At the end, we embed the Laplacian into RG_star, the group ring
# with the same basis as RG but with twisted multiplciation, i.e.
# (1+g)(1+h)=1+g+h+g^(-1)h. This is needed to translate hermitian squares
# to the standard ones for the definition of the semi-definite optimization problem
# (solvers prefer standard squares to hermitian :).
# This has no effect on the shape of the Laplacian since we just embed it.
RG_Δ₃_star = LowCohomologySOS.group_ring(sln, half_basis_Δ[3], star_multiplication = true)
RG_Δ₄_star = LowCohomologySOS.group_ring(sln, half_basis_Δ[4], star_multiplication = true)
RG_Δ₅_star = LowCohomologySOS.group_ring(sln, half_basis_Δ[5], star_multiplication = true)
RG_Δ₆_star = LowCohomologySOS.group_ring(sln, half_basis_Δ[6], star_multiplication = true)
RG_Δ₇_star = LowCohomologySOS.group_ring(sln, half_basis_Δ[7], star_multiplication = true)

d₄x = LowCohomologySOS.embed.(identity, d₄, Ref(RG_Δ[3]));
Δ₃x = d₄x*d₄x'+dim3_stab_part

d₄x = LowCohomologySOS.embed.(identity, d₄, Ref(RG_Δ[4]));
d₅x = LowCohomologySOS.embed.(identity, d₅, Ref(RG_Δ[4]));
Δ₄x = d₄x'*d₄x+d₅x*d₅x'+dim4_stab_part

d₅x = LowCohomologySOS.embed.(identity, d₅, Ref(RG_Δ[5]));
d₆x = LowCohomologySOS.embed.(identity, d₆, Ref(RG_Δ[5]));
Δ₅x = d₅x'*d₅x+d₆x*d₆x'+dim5_stab_part

d₆x = LowCohomologySOS.embed.(identity, d₆, Ref(RG_Δ[6]));
d₇x = LowCohomologySOS.embed.(identity, d₇, Ref(RG_Δ[6]));
Δ₆x = d₆x'*d₆x+d₇x*d₇x'+dim6_stab_part

d₇x = LowCohomologySOS.embed.(identity, d₇, Ref(RG_Δ[7]));
d₈x = LowCohomologySOS.embed.(identity, d₈, Ref(RG_Δ[7]));
Δ₇x = d₇x'*d₇x+d₈x*d₈x'+dim7_stab_part

Δ₃ = LowCohomologySOS.embed.(identity, Δ₃x, Ref(RG_Δ₃_star))
Δ₄ = LowCohomologySOS.embed.(identity, Δ₄x, Ref(RG_Δ₄_star))
Δ₅ = LowCohomologySOS.embed.(identity, Δ₅x, Ref(RG_Δ₅_star))
Δ₆ = LowCohomologySOS.embed.(identity, Δ₆x, Ref(RG_Δ₆_star))
Δ₇ = LowCohomologySOS.embed.(identity, Δ₇x, Ref(RG_Δ₇_star));

# Check if the Laplacians are hermitian (just to be sure we have not spolied something obvious)
@info Δ₃' == Δ₃
@info Δ₄' == Δ₄
@info Δ₅' == Δ₅
@info Δ₆' == Δ₆
@info Δ₇' == Δ₇

# Load the solver ...
# include("/home/mizerka/Desktop/LowCohomologySOS/scripts/optimizers.jl");
# include("/Users/piotrmizerka/Desktop/postdoc_warsaw/code/LowCohomologySOS/scripts/optimizers.jl");
# ... and the optimization package
# using JuMP
