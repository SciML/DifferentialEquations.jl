######
##FEM Stochastic Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(5)
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet)
prob = prob_poisson_noisywave

sol = solve(fem_mesh::FEMmesh,prob::PoissonProblem)#,solver=:CG) #TODO Fix CG and switch back

TEST_PLOT && plot(sol,title=["True Deterministic Solution" "Stochastic Solution"],plot_analytic=true)
#This condition should be true with really high probability
var(sol.u) < 8e-4
