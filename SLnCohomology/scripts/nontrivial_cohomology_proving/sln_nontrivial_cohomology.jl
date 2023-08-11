# These parameters (hard-coded) are subject to appropriate change.
const N = 3
const p = 3

# include("../differentials_computation/sln_laplacians.jl");
# Instead, load from precomputed Laplacian data:
sln_laplacian_data = deserialize(joinpath(@__DIR__, "../differentials_computation/precomputed_laplacians/sl"*string(N)*"_laplacians.sjl"))
Δ = sln_laplacian_data["laplacians"]
differential_degrees = sln_laplacian_data["differential_degrees"]

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

# Compute π(Δₙ) for π, the representation given by permutation representation of SL(n,p)
# and n varying through homology degrees.
πΔ = Dict()
for entry in Δ
    n = entry[1]
    πΔ[n] = vcat(
        [
            hcat([representing_matrix(Δ[n][i,j],p) for j in 1:size(Δ[n])[2]]...) 
            for i in 1:size(Δ[n])[1]
        ]...
    )
end

# For SL(3,Z), modular projection with p = 3 ##############################################################################
for entry in Δ
    n = entry[1]
    @info "rank H^"*string(differential_degrees[end]-n)*" = "*string(size(πΔ[n])[1]-rank(πΔ[n]))
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
# ##########################################################################################
