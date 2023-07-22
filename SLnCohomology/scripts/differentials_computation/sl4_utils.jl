using Pkg
# Pkg.activate(normpath(joinpath(@__DIR__, "./")))
Pkg.activate("/home/mizerka/Desktop/LowCohomologySOS") # you need to change this path for the path of LowCohomologySOS at your computer
# Pkg.activate("/Users/piotrmizerka/Desktop/postdoc_warsaw/code/LowCohomologySOS")

# We employ the package "Groups" of Marek Kaluba (https://github.com/kalmarek/Groups.jl)
# and the package "LowCohomologySOS" of Marek and Piotr (https://github.com/piotrmizerka/LowCohomologySOS)
# Pkg.add("Groups") # once added, no need to try to add it again - this line can be commented out in this case
# Pkg.add(url="https://github.com/piotrmizerka/LowCohomologySOS") # as above
using Groups
using LowCohomologySOS
using JuMP
using Serialization
using LinearAlgebra

include("../../src/helpful_functions.jl")

const N = 4
sln = MatrixGroups.SpecialLinearGroup{N}(Int8)
sln_gens = gens(sln)

# Load the data from Ben
sl4_bound_stab = deserialize(joinpath(@__DIR__, "./sl4_bound_stab.sjl"))

# Compute the stabilisers as subgroups of SL(n,ℤ)
m3_stabs = [[(gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][3]]
m4_stabs = [[(gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][4]]
m5_stabs = [[(gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][5]]
m6_stabs = [[(gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][6]]
m7_stabs = [[(gelt_from_matrix(M,sln),σ) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][7]]

m3_stabsx = [[gelt_from_matrix(M,sln) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][3]]
m4_stabsx = [[gelt_from_matrix(M,sln) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][4]]
m5_stabsx = [[gelt_from_matrix(M,sln) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][5]]
m6_stabsx = [[gelt_from_matrix(M,sln) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][6]]
m7_stabsx = [[gelt_from_matrix(M,sln) for (M,σ) in stab] for stab in sl4_bound_stab["stabilisers"][7]]

# Compute the support (i.e. half_basis) for the group ring over which we'll try to find a spectral gap.
# half_basis is quite minimalistic - we just add to half_basis the coset elements appearing in the differentials.
d_union = Dict(4=>[one(sln)],5=>[one(sln)],6=>[one(sln)],7=>[one(sln)],8=>[one(sln)])
for k in 4:8
    for i in eachindex(sl4_bound_stab["boundaries"][k])
        for j in eachindex(sl4_bound_stab["boundaries"][k][i])
            coset = sl4_bound_stab["boundaries"][k][i][j]["orbit_coset_with_orientation"]
            d_union[k] = union(d_union[k],[gelt_from_matrix(M,sln) for (M,σ) in coset])
        end
    end
end
half_basis_Δ = Dict()
half_basis_Δ[3] = d_union[4]
half_basis_Δ[3] = unique([half_basis_Δ[3];inv.(half_basis_Δ[3])])
for k in 5:8
    half_basis_Δ[k-1] = union(d_union[k-1], d_union[k])
    half_basis_Δ[k-1] = unique([half_basis_Δ[k-1];inv.(half_basis_Δ[k-1])])
end

# Compute the group rings (with standard multiplciation by convolution: (1+g)(1+h)=1+g+h+gh)
RG_Δ = Dict()
RG_Δ[3] = LowCohomologySOS.group_ring(sln, half_basis_Δ[3], star_multiplication = false)
RG_Δ[4] = LowCohomologySOS.group_ring(sln, half_basis_Δ[4], star_multiplication = false)
RG_Δ[5] = LowCohomologySOS.group_ring(sln, half_basis_Δ[5], star_multiplication = false)
RG_Δ[6] = LowCohomologySOS.group_ring(sln, half_basis_Δ[6], star_multiplication = false)
RG_Δ[7] = LowCohomologySOS.group_ring(sln, half_basis_Δ[7], star_multiplication = false)

# Compute the differential as matrices over RGs
cells_number = Dict(3=>1, 4=>3, 5=>4, 6=>4, 7=>2, 8=>2)
d = Dict([(k,i,j) => 0//1*zero(RG_Δ[k-1]) for k in 4:8 for i in 1:cells_number[k-1] for j in 1:cells_number[k]])
for k in 4:8
    for j in eachindex(sl4_bound_stab["boundaries"][k])
        for summand in sl4_bound_stab["boundaries"][k][j]
            coset_orient = [(gelt_from_matrix(M,sln),σ) for (M,σ) in summand["orbit_coset_with_orientation"]]
            summand_in_RG = summand["sign"]*averaged_rep(coset_orient,RG_Δ[k-1])
            d[k,summand["orbit_standard_cell"],j] += summand_in_RG
        end
    end
end
d₄ = [d[4,i,j] for i in 1:cells_number[3], j in 1:cells_number[4]]
d₅ = [d[5,i,j] for i in 1:cells_number[4], j in 1:cells_number[5]]
d₆ = [d[6,i,j] for i in 1:cells_number[5], j in 1:cells_number[6]]
d₇ = [d[7,i,j] for i in 1:cells_number[6], j in 1:cells_number[7]]
d₈ = [d[8,i,j] for i in 1:cells_number[7], j in 1:cells_number[8]];
