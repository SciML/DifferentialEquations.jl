######
##FEM Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
pdeProb = poissonProblemExample_wave()

res = solve(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="CG")

solplot_appxvstrue(res) #To save the plot, add ,savefile="plot.png" or "plot.pdf", etc.

#Error should be low
res.errors["L2"] < 1e-4
