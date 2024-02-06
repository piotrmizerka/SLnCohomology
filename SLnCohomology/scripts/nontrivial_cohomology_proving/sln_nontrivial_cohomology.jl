# These parameters are subject to appropriate change.
const N = 4
const p = 2

# include("../differentials_computation/sln_laplacians.jl");
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using AbstractAlgebra
using BlockArrays
using GAP
using Groups
using LowCohomologySOS
using JSON
using Permutations
using Serialization
using SLnCohomology
using SparseArrays

sln_laplacian_data = deserialize(joinpath(@__DIR__, "../differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
Δ = sln_laplacian_data["laplacians"]
if N == 4 # don't consider 1st cohomology for SL(4,Z) since one of the cells was not simplicial
    delete!(Δ,8)
end

# Compute orthogonal representations without invariant vectors of a subgroup H of SL(N,p).
subgroup_rep = Dict()
if N == 3
    # Linear rep of order 2 of Q₈ embedded in SL(3,3).
    # Elements of Q₈ grouped by conjugacy classes:
    Q8 = [
        [1 0 0;0 1 0;0 0 1], # (order 1),
        [2 0 0;0 2 0;0 0 1], # (order 2),
        [2 2 0;2 1 0;0 0 1], [1 1 0;1 2 0;0 0 1], # (order 4),
        [0 2 0;1 0 0;0 0 1], [0 1 0;2 0 0;0 0 1], # (order 4),
        [1 2 0;2 2 0;0 0 1], [2 1 0;1 1 0;0 0 1] # (order 4)
    ]
    function order_two_q8(g)
        if g == [0 2 0;1 0 0;0 0 1] || g == [0 1 0;2 0 0;0 0 1] || 
           g == [2 2 0;2 1 0;0 0 1] || g == [1 1 0;1 2 0;0 0 1]
            return reshape([-1],1,1)
        else
            return reshape([1],1,1)
        end
    end
    for g in Q8
        subgroup_rep[g] = order_two_q8(g)
    end
    deg = 1
elseif N == 4 || N == 5
    GAP.evalstr("G := SL("*string(N)*","*string(p)*");;")
    GAP.evalstr("conj_cl_sbgps_G := ConjugacyClassesSubgroups(G);;")
    if N == 4 && p == 2
        GAP.evalstr("H_in_G := Representative(conj_cl_sbgps_G[132]);;")
    elseif N == 5 && p == 2
        GAP.evalstr("H_in_G := Representative(conj_cl_sbgps_G[1407]);;")
    end
    GAP.evalstr("H_list := AsList(H_in_G);;")
    GAP.evalstr("iso := IsomorphismPermGroup(H_in_G);;")
    permutations = []
    for i in 1:GAP.evalstr("Order(H_in_G);")
        line = string(GAP.evalstr("Image(iso,H_list["*string(i)*"]);"))
        linex = replace(line, r"\(" => ",[")
        linexx = replace(linex, r"\)" => "]")
        json_perm = JSON.parse("["*linexx[2:end]*"]")
        proper_perm = [Int64.(cycle) for cycle in json_perm]
        push!(permutations,proper_perm)
    end
    # For technical reasons (due to Julia permutation conversions), we must write each permutation including trivial cycles
    deg = SLnCohomology.permutations_degree(permutations)
    function integer_matrix(i::Integer)
        return [GAP.evalstr("Int(H_list["*string(i)*"]["*string(k)*","*string(l)*"]);") for k in 1:N,l in 1:N]
    end
    GAP.evalstr("normal_sbgps_H := NormalSubgroups(H_in_G);")
    if N == 4
        GAP.evalstr("index_two_H := normal_sbgps_H[4];") 
    else
        GAP.evalstr("index_two_H := normal_sbgps_H[5];")
    end
    for i in 1:GAP.evalstr("Order(H_in_G);")
        proper_perm = SLnCohomology.standarize_permutation(permutations[i],deg)
        proper_perm_2 = [x for x in proper_perm]
        perm_stand = Permutation(proper_perm_2)
        even_sign = GAP.evalstr("H_list["*string(i)*"] in index_two_H;")
        if !even_sign
            subgroup_rep[integer_matrix(i)] = sparse(-Int.(Matrix(perm_stand)^(-1)))
        else
            subgroup_rep[integer_matrix(i)] = sparse(Int.(Matrix(perm_stand)^(-1)))
        end
    end
end
H = collect(keys(subgroup_rep))
no_inv_subspace = SLnCohomology.no_inv_subspace(H,subgroup_rep)
@assert length(no_inv_subspace) == deg # this means that we have no invariant vectors

# Induce the repesentation of H to SL(N,p)
file_path = joinpath(@__DIR__, "../../scripts/nontrivial_cohomology_proving/sln_p_matrices/sl"*string(N)*"_"*string(p)*"_matrices.txt")
sl_n_p_matrices = SLnCohomology.read_slnp_matrices(file_path,N)
support = SLnCohomology.laplacians_support(Δ,p)
coset_data = SLnCohomology.coset_data(H,sl_n_p_matrices,p)
π = SLnCohomology.ind_rep_dict(support,subgroup_rep,coset_data,deg,p)

# Compute Laplacians evaluated on the induced representation
sbgp_ind = length(coset_data["cosets_representatives"])
πΔ = Dict()
for entry in Δ
    n = entry[1]
    @info n
    πΔ[n] = vcat(
        [
            hcat([SLnCohomology.representing_matrix(Δ[n][i,j],π, deg, p, sbgp_ind) for j in 1:size(Δ[n])[2]]...)
            for i in 1:size(Δ[n])[1]
        ]...
    )
end

# Show cohomologies' ranks for the above Laplacians
for entry in Δ
    n = entry[1]
    @info "rank H^"*string(div(N*(N-1),2)+N-n-1)*" = "*string(size(πΔ[n])[1]-rank(πΔ[n]))
end

# # TODO: move this function to tests
# function representing_matrix_trivial(ξ,p::Integer)
#     result = zero(typeof(first(ξ.coeffs)))
#     RG = parent(ξ)
#     for x in RG.basis
#         result += ξ(x)
#     end
#     return result
# end
# # TODO: move all sanity checks to tests!
# # A sanity check - since H¹(SL(3,Z),π) = 0, we shall get the full rank for the perm repr of Δ₄:
# @assert πΔ[4]' == πΔ[4]
# @assert size(πΔ[4])[1]-rank(πΔ[4]) == 0

# # A sanity check: since the permutation representation has a nontrivial fixed pont set, 
# # we shall get nontrivial zero cohomology for this representation.
# @assert πΔ[5]' == πΔ[5]
# size(πΔ[5])[1]-rank(πΔ[5]) # non-full rank - means nontrivial 0-cohomology

# # H² - we get nontrivial cohomology!
# @assert πΔ[3]' == πΔ[3]
# size(πΔ[3])[1]-rank(πΔ[3]) # non-full rank - means nontrivial 2-cohomology

# @assert πΔ[2]' == πΔ[2]
# size(πΔ[2])[1]-rank(πΔ[2])
