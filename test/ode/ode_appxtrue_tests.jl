using DifferentialEquations

f = (t,u) -> u
prob = ODEProblem(f,1/2)
analytic = (t,u₀) -> u₀*exp(t)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:Euler)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(10),alg=:Vern9)

prob2 = ODEProblem(f,1/2,analytic=analytic)
sol3 =solve(prob_ode_linear,[0,1];Δt=1//2^(4),alg=:Euler)

appxTrue!(sol,sol2)

sol4 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:Euler)
test_sol = TestSolution(sol2)
appxTrue!(sol4,test_sol)

sol.appxTrue == true && sol.errors[:L2] ≈ 0.018865798306718855 && sol.errors[:L2] ≈ sol4.errors[:L2]
