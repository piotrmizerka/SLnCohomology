# TODO: verify the unified new code with the old one!!!

# Compare two matrices over RSlns which don't necessarily have the same underyling
# group rings. We use matrix evaluation of elts of two unerlying slns for comparison.
using SparseArrays
using StaticArrays
function compare_rg_matrices(A,B)
    if size(A) != size(B)
        return false
    end
    RG_A = parent(first(A))
    RG_B = parent(first(B))
    sln_A = parent(first(RG_A.basis))
    sln_B = parent(first(RG_B.basis))
    for i in 1:size(A)[1]
        for j in 1:size(A)[2]
            aij_nonzero_inds = MArray(SparseArrays.nonzeroinds(A[i,j].coeffs))
            bij_nonzero_inds = MArray(SparseArrays.nonzeroinds(B[i,j].coeffs))
            if length(aij_nonzero_inds) != length(bij_nonzero_inds)
                return false
            end
            for ind in aij_nonzero_inds
                # g_A = RG_A.basis[ind]
                # @info MatrixGroups.matrix_repr(g_A)
                @info SLnCohomology.gelt_from_matrix(MatrixGroups.matrix_repr(RG_A.basis[ind]),sln_A)
                # g_B = SLnCohomology.gelt_from_matrix(MatrixGroups.matrix_repr(g_A),sln_A)
                # g_B = SLnCohomology.gelt_from_matrix(Int8[0 1 0; 0 0 1; 1 1 0],sln_B)
                # if A[i,j](g_A) != B[i,j](g_B)
                #     return false
                # end
            end
        end
    end
    return true
end
compare_rg_matrices(d[3],d₃)
size(d[3])
RG1 = parent(d[3][1,1])
RG2 = parent(d₃[1,1])
parent(RG1.basis[1])
parent(RG2.basis[1])
RG2.basis
d₃[1,1](one(sl3))
d[3][1,1](one(sln))
RG1 = parent(d[5][1,1])
RG2 = parent(d₅[1,1])
parent(RG1.basis[1])
parent(RG2.basis[1])
RG2.basis
for g in RG1.basis
    # @assert d[3][1,1](g) == d₃[1,1](g)
    @assert d[3][1,1](g) == d[3][1,1](g)
end
d[4]
d₄
d[5]
d₅

dim2_stab_part[1,1]
stab_part_dim[2][1,1]
Δ[2][1,1]
Δ₂[1,1](one(sl3))
Δ[3][1,1]
Δ₃[1,1]
Δ[4]
Δ₄
Δ[5]
Δ₅