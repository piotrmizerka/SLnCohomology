# Verify the new code with the old using sl3 data.

# Compare two matrices over RSlns which don't necessarily have the same underyling
# group rings. We use matrix evaluation of elts of two unerlying slns for comparison.
include("sln_laplacians.jl");
include("sl3_laplacians.jl");

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))

using SparseArrays

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
            aij_nonzero_inds = SparseArrays.nonzeroinds(A[i,j].coeffs)
            bij_nonzero_inds = SparseArrays.nonzeroinds(B[i,j].coeffs)
            if length(aij_nonzero_inds) != length(bij_nonzero_inds)
                return false
            end
            for ind in aij_nonzero_inds
                g_A = RG_A.basis[ind]
                g_A_matrix = Array(MatrixGroups.matrix_repr(g_A))
                g_B = SLnCohomology.gelt_from_matrix(g_A_matrix,sln_B)
                if A[i,j](g_A) != B[i,j](g_B)
                    return false
                end
            end
        end
    end
    return true
end

compare_rg_matrices(d[3],d₃)
compare_rg_matrices(d[4],reshape(d₄,2,1))
compare_rg_matrices(d[5],d₅)

compare_rg_matrices(Δ[2],Δ₂)
compare_rg_matrices(Δ[3],Δ₃)
compare_rg_matrices(Δ[4],Δ₄)
compare_rg_matrices(Δ[5],Δ₅)
