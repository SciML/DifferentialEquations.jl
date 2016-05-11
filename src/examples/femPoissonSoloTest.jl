######
##FEM Poisson Method Tests
######
using DiffEq

T = 1
Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
pdeProb = poissonProblemExample_wave()

res = fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="GMRES")
solplot_appxvstrue(res,savefile="plot.svg")
