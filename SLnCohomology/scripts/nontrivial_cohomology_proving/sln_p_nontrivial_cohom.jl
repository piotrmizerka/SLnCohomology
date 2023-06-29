# include("../differentials_computation/sl3_laplacians.jl");
include("../differentials_computation/sl4_laplacians.jl");

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using Permutations
using JSON
using GAP

const p = 2
const N = 4

function slnp_order(n,p)
    result = 1
    for i in 0:(n-1)
        result *= (p^n-p^i)
    end
    result = div(result,p-1)
    return result
end

# Compute a permutation representation for SL(N,p) and store it as a dictionary "permutation_matrices"
# which assigns a given matrix from SL(N,p) a permutation matrix. #############################################################################

# Before running the code below, requires running the script "permutation_matrices.g" in gap which saves the elements of SL(N,p) 
# as matrices in a separate file "./sln_p_matrices.txt" via GAP script "permutation_matrices.g".
# One has to do this from the command line and also one has to change the parameters N and p in "permutation_matrices.g"
# manually (sorry for the hard-coding :)). The following commands do this (run the terminal in the same directory as this julia file):
# gap
# gap> Read(".sln_p_matrices.g");
# Note that some tiny cleaning-up of "./sln_p_matrices.txt" may be needed so that it contains no words like "gap" or symbols ">", etc.
# I did this manually - these lines appeared only at the beginning and at the end apparently.

# REMARK: I pushed the data for SL(3,3) to the online repo, so there is no need to do that for that case. 

# Compute permutations for each matrix from SL(N,p) with GAP (we use a Julia wrapper for GAP)
GAP.evalstr("G := SL("*string(N)*","*string(p)*");;")
GAP.evalstr("G_list := AsList(G);;")
GAP.evalstr("iso := IsomorphismPermGroup(G);;")
permutations = []
for i in 1:slnp_order(N,p)
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
                deg = maximum(cycle)
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

# file_path = "./sln_p_matrices.txt"
file_path = "./sl4_2_matrices.txt"
file = open(file_path, "r")
i = 0
current_matrix = Int8.(zeros(N,N))
permutation_matrices = Dict()
for line in eachline(file)
    linex = replace(line, r"\s+" => "")
    linexx = replace(linex, r"\." => "0")
    for j in eachindex(linexx)
        current_matrix[i%N+1,j] = parse(Int8,linexx[j])
    end
    if i%N == N-1
        proper_perm = standarize_permutation(permutations[div(i,N)+1], deg)
        proper_perm_2 = [x for x in proper_perm]
        permutation_matrices[current_matrix] = Matrix(Permutation(proper_perm_2))
        current_matrix = Int8.(zeros(N,N))
    end
    i += 1
end
close(file)
# rm("./sln_p_matrices.txt") don't erase the elts of SL(N,p) at this point
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
    for x in RG.basis
        proj = projection(MatrixGroups.matrix_repr(x),p)
        result += ξ(x)*permutation_matrices[proj]
    end
    return result
end

# For SL(3,Z), modular projection with p = 3 ##############################################################################
# A sanity check - since H^1(SL(3,Z),π) = 0, we shall get the full rank for the perm repr of Δ₄:
delta_4_perm = representing_matrix(Δ₄[1,1],p)
@assert delta_4_perm' == delta_4_perm
@assert size(delta_4_perm)[1]-rank(delta_4_perm) == 0

# A sanity check: since the permutation representation has a nontrivial fixed pont set, 
# we shall get nontrivial zero cohomology for this representation.
delta_5_perm = representing_matrix(Δ₅[1,1],p)
@assert delta_5_perm' == delta_5_perm
size(delta_5_perm)[1]-rank(delta_5_perm) # non-full rank - means nontrivial 0-cohomology

# H^2 - we get nontrivial cohomology!
delta_3_11 = representing_matrix(Δ₃[1,1],p)
delta_3_12 = representing_matrix(Δ₃[1,2],p)
delta_3_21 = representing_matrix(Δ₃[2,1],p)
delta_3_22 = representing_matrix(Δ₃[2,2],p)
delta_3_perm = [hcat(delta_3_11,delta_3_12);hcat(delta_3_21,delta_3_22)]
@assert delta_3_perm' == delta_3_perm
size(delta_3_perm)[1]-rank(delta_3_perm) # non-full rank - means nontrivial 2-cohomology
##########################################################################################

# For SL(4,Z), modular projection with p = 2 ##############################################################################
# A sanity check - since, by Bader-Sauer, H^2(SL(4,Z),π) = 0,
#  we shall get the full rank for the perm repr of Δ₄:
delta_7_11 = representing_matrix(Δ₇[1,1],p)
delta_7_12 = representing_matrix(Δ₇[1,2],p)
delta_7_21 = representing_matrix(Δ₇[2,1],p)
delta_7_22 = representing_matrix(Δ₇[2,2],p)
delta_7_perm = [hcat(delta_7_11,delta_7_12);hcat(delta_7_21,delta_7_22)]
@assert delta_7_perm' == delta_7_perm
@assert size(delta_7_perm)[1]-rank(delta_7_perm) == 0

# H^3 - we get nontrivial cohomology!
delta_6_11 = representing_matrix(Δ₆[1,1],p)
delta_6_12 = representing_matrix(Δ₆[1,2],p)
delta_6_13 = representing_matrix(Δ₆[1,3],p)
delta_6_14 = representing_matrix(Δ₆[1,4],p)
delta_6_21 = representing_matrix(Δ₆[2,1],p)
delta_6_22 = representing_matrix(Δ₆[2,2],p)
delta_6_23 = representing_matrix(Δ₆[2,3],p)
delta_6_24 = representing_matrix(Δ₆[2,4],p)
delta_6_31 = representing_matrix(Δ₆[3,1],p)
delta_6_32 = representing_matrix(Δ₆[3,2],p)
delta_6_33 = representing_matrix(Δ₆[3,3],p)
delta_6_34 = representing_matrix(Δ₆[3,4],p)
delta_6_41 = representing_matrix(Δ₆[4,1],p)
delta_6_42 = representing_matrix(Δ₆[4,2],p)
delta_6_43 = representing_matrix(Δ₆[4,3],p)
delta_6_44 = representing_matrix(Δ₆[4,4],p)
delta_6_perm = [
    hcat(delta_6_11,delta_6_12,delta_6_13,delta_6_14);
    hcat(delta_6_21,delta_6_22,delta_6_23,delta_6_24);
    hcat(delta_6_31,delta_6_32,delta_6_33,delta_6_34);
    hcat(delta_6_41,delta_6_42,delta_6_43,delta_6_44)
]
@assert delta_6_perm' == delta_6_perm
size(delta_6_perm)[1]-rank(delta_6_perm)

# H^4 - here our perm rep yields trivial cohomology
delta_5_11 = representing_matrix(Δ₅[1,1],p)
delta_5_12 = representing_matrix(Δ₅[1,2],p)
delta_5_13 = representing_matrix(Δ₅[1,3],p)
delta_5_14 = representing_matrix(Δ₅[1,4],p)
delta_5_21 = representing_matrix(Δ₅[2,1],p)
delta_5_22 = representing_matrix(Δ₅[2,2],p)
delta_5_23 = representing_matrix(Δ₅[2,3],p)
delta_5_24 = representing_matrix(Δ₅[2,4],p)
delta_5_31 = representing_matrix(Δ₅[3,1],p)
delta_5_32 = representing_matrix(Δ₅[3,2],p)
delta_5_33 = representing_matrix(Δ₅[3,3],p)
delta_5_34 = representing_matrix(Δ₅[3,4],p)
delta_5_41 = representing_matrix(Δ₅[4,1],p)
delta_5_42 = representing_matrix(Δ₅[4,2],p)
delta_5_43 = representing_matrix(Δ₅[4,3],p)
delta_5_44 = representing_matrix(Δ₅[4,4],p)
delta_5_perm = [
    hcat(delta_5_11,delta_5_12,delta_5_13,delta_5_14);
    hcat(delta_5_21,delta_5_22,delta_5_23,delta_5_24);
    hcat(delta_5_31,delta_5_32,delta_5_33,delta_5_34);
    hcat(delta_5_41,delta_5_42,delta_5_43,delta_5_44)
]
@assert delta_5_perm' == delta_5_perm
size(delta_5_perm)[1]-rank(delta_5_perm)

# H^5:
delta_4_11 = representing_matrixx(Δ₄[1,1],p)
delta_4_12 = representing_matrixx(Δ₄[1,2],p)
delta_4_13 = representing_matrixx(Δ₄[1,3],p)
delta_4_21 = representing_matrixx(Δ₄[2,1],p)
delta_4_22 = representing_matrixx(Δ₄[2,2],p)
delta_4_23 = representing_matrixx(Δ₄[2,3],p)
delta_4_31 = representing_matrixx(Δ₄[3,1],p)
delta_4_32 = representing_matrixx(Δ₄[3,2],p)
delta_4_33 = representing_matrixx(Δ₄[3,3],p)
delta_4_perm = [
    hcat(delta_4_11,delta_4_12,delta_4_13);
    hcat(delta_4_21,delta_4_22,delta_4_23);
    hcat(delta_4_31,delta_4_32,delta_4_33)
]
@assert delta_4_perm' == delta_4_perm
size(delta_4_perm)[1]-rank(delta_4_perm)
##########################################################################################