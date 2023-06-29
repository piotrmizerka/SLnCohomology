include("../differentials_computation/sl4_laplacians.jl");

# Δ₅
I = [
    one(RG_Δ₅_star) zero(RG_Δ₅_star) zero(RG_Δ₅_star) zero(RG_Δ₅_star);
    zero(RG_Δ₅_star) one(RG_Δ₅_star) zero(RG_Δ₅_star) zero(RG_Δ₅_star);
    zero(RG_Δ₅_star) zero(RG_Δ₅_star) one(RG_Δ₅_star) zero(RG_Δ₅_star);
    zero(RG_Δ₅_star) zero(RG_Δ₅_star) zero(RG_Δ₅_star) one(RG_Δ₅_star)
]
sos_problem = LowCohomologySOS.sos_problem(Δ₅, I, 0.05)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 20_000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem)
LowCohomologySOS.certify_sos_decomposition(Δ₅, I, λ, Q, half_basis_Δ[5])