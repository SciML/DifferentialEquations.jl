######
##FEM Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(3)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Neumann")
pdeProb = poissonProblemExample_birthdeath()

sol = solve(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="GMRES")

plot(sol,plottrue=false,zlim=(0,2))

#Returns true if computed solution is homogenous near 2
maximum(abs(sol.u - 2))< 1e-12
