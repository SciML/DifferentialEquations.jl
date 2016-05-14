######
##FEM Stochastic Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
pdeProb = poissonProblemExample_noisyWave()

res = fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="GMRES")

solplot_appxvstrue(res,trueTitle="True Deterministic Solution",appxTitle="Stochastic Solution")

#This condition should be true with really high probability
var(res.u) < 2e-4
