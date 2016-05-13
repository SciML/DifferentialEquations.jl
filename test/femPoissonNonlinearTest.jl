######
##FEM Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(6)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Neumann")
pdeProb = poissonProblemExample_birthdeath()

res = fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="GMRES")

if !isdefined(:testState) #Don't plot during test
  Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)
  #solplot_appx(res,savefile="plot.svg") #Approximately Flat solution at 2
end
