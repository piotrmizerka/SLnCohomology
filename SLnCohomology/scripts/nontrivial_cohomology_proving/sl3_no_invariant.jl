const N = 3
const p = 3

# include("../differentials_computation/sln_laplacians.jl");
# Instead, load from precomputed Laplacian data: ##########################
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using BlockArrays
using Groups
using LowCohomologySOS
using Serialization
using SLnCohomology

sln_laplacian_data = deserialize(joinpath(@__DIR__, "../differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
Δ = sln_laplacian_data["laplacians"]

using Permutations
using JSON
using GAP
using SparseArrays

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

# Load SL(N,p) matrices
sl_n_p_matrices = []
dir_path = "/Users/piotrmizerka/Desktop/postdoc_warsaw/articles/HigherTSL3/SLnCohomology/" # need to be changed accordingly!
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
H = Q8
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
    result = Int8.(spzeros(length(cosets_representatives),length(cosets_representatives)))
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

# TODO: move the code below to tests
it = 0
f = induced_q8 # toggle if necesssary
for g in sl_n_p_matrices
    @info it
    @assert Matrix(f(g))^(-1) == f(g)' # orthogonality test
    @assert Matrix(f(g))^(-1) == f(inverse_elt(g,sl_n_p_matrices,p))
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
ind = induced_q8
it = 0
for g in sl_n_p_matrices
    represention_dict[g] = ind(g)
    it += 1
end

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

function representing_matrix_induction(ξ,p::Integer) # TODO: don't repeat the code!
    result = typeof(first(ξ.coeffs)).(spzeros(length(cosets_representatives),length(cosets_representatives)))
    RG = parent(ξ)
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ(x)*represention_dict[proj]
    end
    return Matrix(result)
end

# Compute π(Δₙ) for n varying through homology degrees.
πΔ = Dict()
for entry in Δ
    n = entry[1]
    @info n
    πΔ[n] = vcat(
        [
            hcat([representing_matrix_induction(Δ[n][i,j],p) for j in 1:size(Δ[n])[2]]...)
            for i in 1:size(Δ[n])[1]
        ]...
    )
end

for entry in Δ
    n = entry[1]
    @info "rank H^"*string(div(N*(N-1),2)+N-n-1)*" = "*string(size(πΔ[n])[1]-rank(πΔ[n]))
end