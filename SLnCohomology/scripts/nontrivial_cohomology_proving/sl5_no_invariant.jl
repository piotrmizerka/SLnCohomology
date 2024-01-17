# TODO: remove reapeated code

const N = 5
const p = 2

# using Pkg
# Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using AbstractAlgebra
using BlockArrays
using GAP
using Groups
using JSON
using LowCohomologySOS
using Permutations
using Serialization
using SLnCohomology
using SparseArrays

sln_laplacian_data = deserialize(joinpath(@__DIR__, "scripts/differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
Δ = sln_laplacian_data["laplacians"]

# TODO: pay attention if no deletes are necessary due to non-standard differentials

# Compute flip-permutation representation of the subgroup of SL(5,2)
# isomorphic to Smallgroup(576,8654), embedded as subgroup of SL(4,2) (cannonically embedded in SL(5,2)). 
# Then we'll induce our π from it. #################################################
GAP.evalstr("G := SL("*string(N)*","*string(p)*");;")
GAP.evalstr("conj_cl_sbgps_sl_5_2 := ConjugacyClassesSubgroups(G);;")
GAP.evalstr("H_in_G := Representative(conj_cl_sbgps_sl_5_2[1407]);;")
GAP.evalstr("H_list := AsList(H_in_G);;")
GAP.evalstr("iso := IsomorphismPermGroup(H_in_G);;")
permutations = []
for i in 1:64512
    line = string(GAP.evalstr("Image(iso,H_list["*string(i)*"]);"))
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

function integer_matrix(i::Integer)
    n = 5
    return [GAP.evalstr("Int(H_list["*string(i)*"]["*string(k)*","*string(l)*"]);") for k in 1:n,l in 1:n]
end

GAP.evalstr("normal_sbgps_H := NormalSubgroups(H_in_G);")
GAP.evalstr("index_two_H := normal_sbgps_H[5];") # this is this group: https://people.maths.bris.ac.uk/~matyd/GroupNames/288i1/PSO+(4,3).html
flip_permutation_matrices = Dict()
for i in 1:64512
    proper_perm = standarize_permutation(permutations[i],deg)
    proper_perm_2 = [x for x in proper_perm]
    perm_stand = Permutation(proper_perm_2)
    even_sign = GAP.evalstr("H_list["*string(i)*"] in index_two_H;")
    if !even_sign
        flip_permutation_matrices[integer_matrix(i)] = sparse(-Int.(Matrix(perm_stand)^(-1)))
    else
        flip_permutation_matrices[integer_matrix(i)] = sparse(Int.(Matrix(perm_stand)^(-1)))
    end
end

sl_n_p_matrices = []
dir_path = "/users/piotrmizerka/HigherTSL3/SLnCohomology/" # need to be changed accordingly!
# dir_path = "/Users/piotrmizerka/Desktop/postdoc_warsaw/articles/HigherTSL3/SLnCohomology/"
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

function matrix_mod_p(M,p::Integer)
    result = copy(M)
    for i in eachindex(result)
        result[i] %= p
    end
    return result
end

cosets = Dict()
cosets_representatives = []
considered_matrices = Dict()
for g in sl_n_p_matrices
    considered_matrices[g] = false
end
H = collect(keys(flip_permutation_matrices))

# Check if what we have has no invariant vectors:
function no_inv(M)
    result = []
    for i in 1:size(M)[1]
        if M[i,i] == -1
            push!(result,i)
        end
    end
    return result
end
no_inv_subspace = Set([])
for h in H
    global no_inv_subspace = union!(no_inv_subspace,no_inv(flip_permutation_matrices[h]))
end
@assert length(no_inv_subspace) == deg # this means that we have no invariant vectors

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

function induced_H(g)
    block_sizes = [deg for i in 1:length(cosets_representatives)]
    total_size = deg*length(cosets_representatives)
    result = BlockArray{Int8}(spzeros(total_size,total_size), block_sizes, block_sizes)
    for i in eachindex(cosets_representatives)
        gi = cosets_representatives[i]
        ggi = matrix_mod_p(g*gi,p)
        gj = cosets[ggi]
        j = cosets_representatives_indices[gj]
        gj_inv = inv(AbstractAlgebra.matrix(GF(p),gj))
        h = Int8.(lift.(gj_inv*AbstractAlgebra.matrix(GF(p),ggi)))
        setblock!(result, flip_permutation_matrices[h], j,i)
    end
    return result
end

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

support = []
for entry in Δ
    n = entry[1]
    @info n
    RG = parent(first(Δ[n]))
    for k in eachindex(Δ[n])
        for i in SparseArrays.nonzeroinds(Δ[n][k].coeffs)
            x = RG.basis[i]
            proj = projection(MatrixGroups.matrix_repr(x),p)
            push!(support,proj)
        end
    end
end
support = collect(Set(support))

# Works up to this point. The following loop takes too much time...
π = Dict()
it = 0
for g in support
    if it == 100
        @info it, "pi dict"
    end
    π[g] = induced_H(g)
    @info it
    it += 1
end

function representing_matrix(ξ,p::Integer) # TODO: don't repeat the code!
    n = deg*length(cosets_representatives)
    result = typeof(first(ξ.coeffs)).(spzeros(n,n))
    RG = parent(ξ)
    for i in SparseArrays.nonzeroinds(ξ.coeffs)
        x = RG.basis[i]
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ(x)*π[proj]
    end
    return Matrix(result)
end

# Compute π(Δₙ) for π, the representation given by permutation representation of SL(n,p)
# and n varying through homology degrees.
πΔ = Dict()
for entry in Δ
    n = entry[1]
    @info n
    πΔ[n] = vcat(
        [
            hcat([representing_matrix(Δ[n][i,j],p) for j in 1:size(Δ[n])[2]]...)
            for i in 1:size(Δ[n])[1]
        ]...
    )
end
# ########################################################################################

for entry in Δ
    n = entry[1]
    @info "rank H^"*string(div(N*(N-1),2)+N-n-1)*" = "*string(size(πΔ[n])[1]-rank(πΔ[n]))
end