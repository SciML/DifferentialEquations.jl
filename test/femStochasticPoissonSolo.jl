######
##FEM Stochastic Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
prob = poissonProblemExample_noisyWave()

sol = solve(femMesh::FEMmesh,prob::PoissonProblem,solver="CG")

plot(sol,title=["True Deterministic Solution" "Stochastic Solution"],plottrue=true)
#This condition should be true with really high probability
var(sol.u) < 8e-4
