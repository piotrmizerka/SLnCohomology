include("../differentials_computation/sl4_laplacians.jl");

# Δ₄
I = [
    one(RG_Δ₄_star) zero(RG_Δ₄_star) zero(RG_Δ₄_star);
    zero(RG_Δ₄_star) one(RG_Δ₄_star) zero(RG_Δ₄_star);
    zero(RG_Δ₄_star) zero(RG_Δ₄_star) one(RG_Δ₄_star)
]
sos_problem = LowCohomologySOS.sos_problem(Δ₄, I, 0.05)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 20_000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem)
LowCohomologySOS.certify_sos_decomposition(Δ₄, I, λ, Q, half_basis_Δ[4])