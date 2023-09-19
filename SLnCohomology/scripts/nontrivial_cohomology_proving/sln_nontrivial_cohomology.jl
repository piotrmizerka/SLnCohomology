# These parameters (hard-coded) are subject to appropriate change.
const N = 4
const p = 2

# include("../differentials_computation/sln_laplacians.jl");
# Instead, load from precomputed Laplacian data: ##########################
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using LowCohomologySOS
using Serialization
using SLnCohomology

sln_laplacian_data = deserialize(joinpath(@__DIR__, "../differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
Δ = sln_laplacian_data["laplacians"]
##########################################################################

using Permutations
using JSON
using GAP
using SparseArrays

# Compute a permutation representation for SL(N,p) and store it as a dictionary "permutation_matrices"
# which assigns a given matrix from SL(N,p) a permutation matrix. #############################################################################

# Before running the code below, requires running the script "sln_p_matrices.g" in gap which saves the elements of SL(N,p) 
# as matrices in a separate file "./slN_p_matrices.txt" (N and p as loaded earlier).
# One has to do this from the command line and also one has to change the parameters N and p in "sln_p_matrices.g"
# manually (sorry for the hard-coding :)). The following commands do this (run the terminal in the same directory as this julia file):
# gap
# gap> Read(".sln_p_matrices.g");
# Note that some tiny cleaning-up of "./slN_p_matrices.txt" may be needed so that it contains no words like "gap" or symbols ">", etc.
# I did this manually - these lines appeared only at the beginning and at the end apparently.

# REMARK: I pushed the data for SL(3,3), SL(4,2), and SL(5,2) to the online repo, 
# so there is no need to do that for these cases. 

# Compute permutations for each matrix from SL(N,p) with GAP (we use a Julia wrapper for GAP)
GAP.evalstr("G := SL("*string(N)*","*string(p)*");;")
GAP.evalstr("G_list := AsList(G);;")
GAP.evalstr("iso := IsomorphismPermGroup(G);;")
permutations = []
for i in 1:SLnCohomology.slnp_order(N,p)
    line = string(GAP.evalstr("Image(iso,G_list["*string(i)*"]);"))
    linex = replace(line, r"\(" => ",[")
    linexx = replace(linex, r"\)" => "]")
    json_perm = JSON.parse("["*linexx[2:end]*"]")
    proper_perm = [Int64.(cycle) for cycle in json_perm]
    push!(permutations,proper_perm)
end

# For technical reasons (due to Julia permutation conversions), we must write each permutation including trivial cycles
deg = 0
for perm in permutations
    for cycle in perm
        if length(cycle) > 0
            if deg < maximum(cycle)
                global deg = maximum(cycle)
            end
        end
    end
end
function standarize_permutation(perm, deg)
    result = (perm[1] == [] ? [] : perm)
    considered = [false for i in 1:deg]
    for cycle in perm
        for x in cycle
            considered[x] = true
        end
    end
    missing_numbers = []
    for i in 1:deg
        if !considered[i]
            push!(missing_numbers,i)
        end
    end
    for x in missing_numbers
        push!(result,[x])
    end
    return result
end

# dir_path = "/home/mizerka/Desktop/HigherTSL3/SLnCohomology" # need to be changed accordingly!
dir_path = "/Users/piotrmizerka/Desktop/postdoc_warsaw/articles/HigherTSL3/SLnCohomology/"
file_path = dir_path*"/scripts/nontrivial_cohomology_proving/sln_p_matrices/sl"*string(N)*"_"*string(p)*"_matrices.txt"
file = open(file_path, "r")
i = 0
current_matrix = Int8.(zeros(N,N))
permutation_matrices = Dict()
for line in eachline(file)
    linex = replace(line, r"\s+" => "")
    linexx = replace(linex, r"\." => "0")
    for j in eachindex(linexx)
        global current_matrix[i%N+1,j] = parse(Int8,linexx[j])
    end
    if i%N == N-1
        proper_perm = standarize_permutation(permutations[div(i,N)+1], deg)
        proper_perm_2 = [x for x in proper_perm]
        permutation_matrices[current_matrix] = Matrix(Permutation(proper_perm_2))
        global current_matrix = Int8.(zeros(N,N))
    end
    global i += 1
end
close(file)
#######################################################################################################################################

# Compute projection onto SL(N,p) from a matrix from SL(N,Z)
function projection(M, p::Integer)
    result = Int8.(zeros(size(M)[1],size(M)[2]))
    for i in 1:size(M)[1]
        for j in 1:size(M)[2]
            result[i,j] = Int8(M[i,j]%p)
            result[i,j] = (result[i,j] >= 0 ? result[i,j] : result[i,j]+p)
        end
    end
    return result
end

# Permutation matrix corresponding do a group ring elt ξ given by the perm repr
function representing_matrix(ξ,p::Integer)
    result = typeof(first(ξ.coeffs)).(zeros(deg,deg))
    RG = parent(ξ)
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ(x)*permutation_matrices[proj]
    end
    return result
end

# TODO: move this function to tests
function representing_matrix_trivial(ξ,p::Integer)
    result = zero(typeof(first(ξ.coeffs)))
    RG = parent(ξ)
    for x in RG.basis
        result += ξ(x)
    end
    return result
end

# Represention without nontrivial invariant vectors - induced from rotations by 180 degrees
function matrix_mod_p(M,p::Integer)
    result = copy(M)
    for i in eachindex(result)
        result[i] %= p
    end
    return result
end
function inverse_elt(g,G,p::Integer)
    I = [i==j ? 1 : 0 for i in 1:size(g)[1],j in 1:size(g)[2]]
    for h in G
        if matrix_mod_p(g*h,p) == I
            return h
        end
    end
end
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
        return -1 #reshape([-1],1,1)
    else
        return 1 #reshape([1],1,1)
    end
end

# TODO: move the code below to tests
for g in Q8
    @assert order_two_q8(g)^(-1) == order_two_q8(inverse_elt(g,Q8,3))
    for h in Q8
        gh = matrix_mod_p(g*h,3)
        @assert order_two_q8(g)*order_two_q8(h) == order_two_q8(gh)
    end
end

# Linear rep of order 2 of S₄ embedded in SL(4,2).
# Elements of S₄ grouped by conjugacy classes:
S4 = [
    # order 1:
    Matrix(Permutation([[1],[2],[3],[4]])),
    # order 2:
    Matrix(Permutation([[1,2],[3],[4]])), Matrix(Permutation([[1,3],[2],[4]])), Matrix(Permutation([[1,4],[2],[3]])), 
    Matrix(Permutation([[2,3],[1],[4]])), Matrix(Permutation([[2,4],[1],[3]])), Matrix(Permutation([[3,4],[1],[2]])),
    # order 2:
    Matrix(Permutation([[1,2],[3,4]])), Matrix(Permutation([[1,3],[2,4]])), Matrix(Permutation([[1,4],[2,3]])),
    # order 3:
    Matrix(Permutation([[1,2,3],[4]])), Matrix(Permutation([[1,3,2],[4]])), 
    Matrix(Permutation([[2,3,4],[1]])), Matrix(Permutation([[2,4,3],[1]])), 
    Matrix(Permutation([[3,4,1],[2]])), Matrix(Permutation([[3,1,4],[2]])),
    Matrix(Permutation([[4,1,2],[3]])), Matrix(Permutation([[4,2,1],[3]])),
    # order 4:
    Matrix(Permutation([[1,2,3,4]])), Matrix(Permutation([[1,2,4,3]])),
    Matrix(Permutation([[1,3,2,4]])), Matrix(Permutation([[1,3,4,2]])),
    Matrix(Permutation([[1,4,2,3]])), Matrix(Permutation([[1,4,3,2]]))
]
function order_two_s4(g)
    if g == Matrix(Permutation([[1,2],[3],[4]])) || g == Matrix(Permutation([[1,3],[2],[4]])) ||
       g == Matrix(Permutation([[1,4],[2],[3]])) || g == Matrix(Permutation([[2,3],[1],[4]])) ||
       g == Matrix(Permutation([[2,4],[1],[3]])) || g == Matrix(Permutation([[3,4],[1],[2]])) ||
       g == Matrix(Permutation([[1,2,3,4]])) || g == Matrix(Permutation([[1,2,4,3]])) ||
       g == Matrix(Permutation([[1,3,2,4]])) || g == Matrix(Permutation([[1,3,4,2]])) ||
       g == Matrix(Permutation([[1,4,2,3]])) || g == Matrix(Permutation([[1,4,3,2]]))
        return -1 #reshape([-1],1,1)
    else
        return 1 #reshape([1],1,1)
    end
end

# TODO: move the code below to tests
for g in S4
    @assert order_two_s4(g)^(-1) == order_two_s4(g^(-1))
    for h in S4
        @assert order_two_s4(g)*order_two_s4(h) == order_two_s4(g*h)
    end
end

# Load SL(N,p) matrices
sl_n_p_matrices = []
# dir_path = "/home/mizerka/Desktop/HigherTSL3/SLnCohomology" # need to be changed accordingly!
dir_path = "/Users/piotrmizerka/Desktop/postdoc_warsaw/articles/HigherTSL3/SLnCohomology/"
file_path = dir_path*"/scripts/nontrivial_cohomology_proving/sln_p_matrices/sl"*string(N)*"_"*string(p)*"_matrices.txt"
file = open(file_path, "r")
i = 0
current_matrix = Int8.(zeros(N,N))
for line in eachline(file)
    linex = replace(line, r"\s+" => "")
    linexx = replace(linex, r"\." => "0")
    for j in eachindex(linexx)
        global current_matrix[i%N+1,j] = parse(Int8,linexx[j])
    end
    if i%N == N-1
        push!(sl_n_p_matrices,current_matrix)
        global current_matrix = Int8.(zeros(N,N))
    end
    global i += 1
end
close(file)

cosets = Dict()
cosets_representatives = []
considered_matrices = Dict()
for g in sl_n_p_matrices
    considered_matrices[g] = false
end
H = S4 # here change for S4 if necesssary
for g in sl_n_p_matrices
    if !considered_matrices[g]
        for h in H
            gh = matrix_mod_p(g*h,p)
            cosets[gh] = g
            considered_matrices[gh] = true
        end
        push!(cosets_representatives,g)
    end
end
cosets_representatives_indices = Dict()
for i in eachindex(cosets_representatives)
    cosets_representatives_indices[cosets_representatives[i]] = i
end

function induced_q8(g)
    result = Int8.(zeros(length(cosets_representatives),length(cosets_representatives)))
    for i in eachindex(cosets_representatives)
        gi = cosets_representatives[i]
        ggi = matrix_mod_p(g*gi,p)
        gj = cosets[ggi]
        j = cosets_representatives_indices[gj]
        for h in H
            if ggi == matrix_mod_p(gj*h,p)
                result[j,i] = order_two_q8(h)
                break
            end
        end
    end
    return result
end

function induced_s4(g)
    result = Int8.(zeros(length(cosets_representatives),length(cosets_representatives)))
    for i in eachindex(cosets_representatives)
        gi = cosets_representatives[i]
        ggi = matrix_mod_p(g*gi,p)
        gj = cosets[ggi]
        j = cosets_representatives_indices[gj]
        for h in H
            if ggi == matrix_mod_p(gj*h,p)
                result[j,i] = order_two_s4(h)
                break
            end
        end
    end
    return result
end

# TODO: move the code below to tests
it = 0
f = induced_s4 # toggle if necesssary
for g in sl_n_p_matrices
    @info it
    @assert f(g)^(-1) == f(inverse_elt(g,sl_n_p_matrices,p))
    it2 = 0
    for h in sl_n_p_matrices
        gh = matrix_mod_p(g*h,p)
        @assert f(g)*f(h) == f(gh)
        if it2 == 10
            break
        end
        it2 +=1 
    end
    if it == 10
        break
    end
    it += 1
end

represention_dict = Dict()
ind = induced_s4
it = 0
for g in sl_n_p_matrices
    @info it
    represention_dict[g] = ind(g)
    it += 1
end

function representing_matrix_induction(ξ,p::Integer) # TODO: don't repeat the code!
    result = typeof(first(ξ.coeffs)).(zeros(length(cosets_representatives),length(cosets_representatives)))
    RG = parent(ξ)
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ(x)*represention_dict[proj] # TODO: sparsify this operation
    end
    return result
end

# Compute π(Δₙ) for π, the representation given by permutation representation of SL(n,p)
# and n varying through homology degrees.
if N == 4 # don't consider 1st cohomology for SL(4,Z) since one of the cells was not simplicial
    # and the orientations were not included for this cell, I guess!
    #  (see the e-mail from Benjamin to all from 2023.03.16)
    delete!(Δ,8)
end

# Δ
# Δ[3]
# Δ[4]
# Δ[5]
# Δ[6]
# Δ[7]
# SparseArrays.nonzeroinds(Δ[3][1,1].coeffs)
# SparseArrays.nonzeroinds(Δ[4][1,1].coeffs)
# SparseArrays.nonzeroinds(Δ[4][1,2].coeffs)
# SparseArrays.nonzeroinds(Δ[4][1,3].coeffs)
# SparseArrays.nonzeroinds(Δ[4][2,1].coeffs)
# SparseArrays.nonzeroinds(Δ[4][2,2].coeffs)
# SparseArrays.nonzeroinds(Δ[4][2,3].coeffs)
# SparseArrays.nonzeroinds(Δ[4][3,1].coeffs)
# SparseArrays.nonzeroinds(Δ[4][3,2].coeffs)
# SparseArrays.nonzeroinds(Δ[4][3,3].coeffs)
# SparseArrays.nonzeroinds(Δ[5][1,1].coeffs)
# SparseArrays.nonzeroinds(Δ[5][2,2].coeffs)
# SparseArrays.nonzeroinds(Δ[5][3,3].coeffs)
# SparseArrays.nonzeroinds(Δ[5][4,4].coeffs)
# SparseArrays.nonzeroinds(Δ[6][1,1].coeffs)
# SparseArrays.nonzeroinds(Δ[6][2,2].coeffs)
# SparseArrays.nonzeroinds(Δ[6][3,3].coeffs)
# SparseArrays.nonzeroinds(Δ[6][4,4].coeffs)
# SparseArrays.nonzeroinds(Δ[7][1,1].coeffs)
# SparseArrays.nonzeroinds(Δ[7][1,2].coeffs)
# SparseArrays.nonzeroinds(Δ[7][2,1].coeffs)
# SparseArrays.nonzeroinds(Δ[7][2,2].coeffs)

πΔ = Dict()
for entry in Δ
    n = entry[1]
    @info n
    πΔ[n] = vcat(
        [
            hcat([representing_matrix_induction(Δ[n][i,j],p) for j in 1:size(Δ[n])[2]]...) # TODO: hard-coded parameter - change!
            for i in 1:size(Δ[n])[1]
        ]...
    )
end

for entry in Δ
    n = entry[1]
    @info "rank H^"*string(div(N*(N-1),2)+N-n-1)*" = "*string(size(πΔ[n])[1]-rank(πΔ[n]))
end

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
