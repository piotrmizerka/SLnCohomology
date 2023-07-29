const N = 3 # SL(N,Z); change accordingly
const n = 2 # Δₙ; change accordingly

include("../differentials_computation/sln_laplacians.jl");

I = [
    i==j ? one(RG_Δ_star[n]) : zero(RG_Δ_star[n]) 
    for i in 1:size(Δ[n])[1], j in 1:size(Δ[n])[2]
]
sos_problem = LowCohomologySOS.sos_problem(Δ[n], I, 0.05) # on the cost of optimality, bound the gap frome above for a quicker solution
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 20_000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem)
LowCohomologySOS.certify_sos_decomposition(Δ₃, I, λ, Q, half_basis_Δ[3])