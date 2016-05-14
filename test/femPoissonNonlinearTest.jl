######
##FEM Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(3)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Neumann")
pdeProb = poissonProblemExample_birthdeath()

res = fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="GMRES")

Plots.gr()
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)

#Returns true if computed solution is homogenous near 2
maximum(abs(res.u - 2))< 1e-12
