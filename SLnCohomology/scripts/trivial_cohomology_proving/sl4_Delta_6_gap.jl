include("../differentials_computation/sl4_laplacians.jl");

# Δ₆
I = [
    one(RG_Δ₆_star) zero(RG_Δ₆_star) zero(RG_Δ₆_star) zero(RG_Δ₆_star);
    zero(RG_Δ₆_star) one(RG_Δ₆_star) zero(RG_Δ₆_star) zero(RG_Δ₆_star);
    zero(RG_Δ₆_star) zero(RG_Δ₆_star) one(RG_Δ₆_star) zero(RG_Δ₆_star);
    zero(RG_Δ₆_star) zero(RG_Δ₆_star) zero(RG_Δ₆_star) one(RG_Δ₆_star)
]
sos_problem = LowCohomologySOS.sos_problem(Δ₆, I, 0.05)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 20_000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem)
LowCohomologySOS.certify_sos_decomposition(Δ₆, I, λ, Q, half_basis_Δ[6])