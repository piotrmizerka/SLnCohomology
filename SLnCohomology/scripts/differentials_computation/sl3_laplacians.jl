# Let's load the precompuded differential data
include("sl3_utils.jl");

# Sanity checks for vanishing compositions of differentials
@assert d₄*d₅ == reshape([zero(RG); zero(RG)],2,1)
@assert d₃*d₄ == [zero(RG)]

# The stabiliser parts which we have to add to get free modules.
# Hopefully (and quite suprisingly for me) the stabilisers' elements belong to half_basis.
dim2_stab_part = reshape([one(RG)-averaged_rep(m2_stab, half_basis, RG)], 1, 1)
dim3_stab_part = [
    one(RG)-averaged_rep(m31_stab, half_basis, RG) zero(RG);
    zero(RG) one(RG)-averaged_rep(m32_stab, half_basis, RG)
]
dim4_stab_part = reshape([one(RG)-averaged_rep(m4_stab, half_basis, RG)], 1, 1)
dim5_stab_part = reshape([one(RG)-averaged_rep(m5_stab, half_basis, RG)], 1, 1);

# Compute the Laplacians.
Δ₂ = d₃*d₃'+dim2_stab_part
Δ₃ = d₃'*d₃+d₄*d₄'+dim3_stab_part
Δ₄ = reshape([d₄'*d₄],1,1)+d₅*d₅'+dim4_stab_part
Δ₅ = d₅'*d₅+dim5_stab_part

# Check if the Laplacians are hermitian (just to be sure we have not spolied something obvious)
@assert Δ₂' == Δ₂
@assert Δ₃' == Δ₃
@assert Δ₄' == Δ₄
@assert Δ₅' == Δ₅

# A sanity check that for the trivial rep H^2 = 0
Δ₃_triv = [0 for i in 1:size(Δ₃)[1], j in 1:size(Δ₃)[2]]
for i in 1:size(Δ₃)[1]
    for j in 1:size(Δ₃)[2]
        Δ₃_triv[i,j] = sum(Δ₃[i, j].coeffs)
    end
end
@assert Δ₃_triv == [1 0;0 1] 
