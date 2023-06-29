include("../differentials_computation/sl4_laplacians.jl");

# Δ₃
I = reshape([one(RG_Δ₃_star)],1,1)
sos_problem = LowCohomologySOS.sos_problem(Δ₃, I, 0.05) # on the cost of optimality, bound the gap frome above for a quicker solution
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 20_000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem)
LowCohomologySOS.certify_sos_decomposition(Δ₃, I, λ, Q, half_basis_Δ[3])