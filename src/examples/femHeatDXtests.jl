######
##FEM Heat Δx Convergence Tests
######

T = 1
Δt = 1//2^(14)
N = 4
topΔx = 7
pdeProb = heatProblemExample_moving()

alg = "Euler"
solutions = cell(N)
for i = 1:N
  Δx = 1//2^(topΔx-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg)
  solutions[i] = res
end
simres = ConvergenceSimulation(solutions)

#solplot_appxvstrue(solutions[1],savefile="plot.svg")

alg = "ImplicitEuler"
solutions = cell(N)
for i = 1:N
  Δx = 1//2^(topΔx-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg)
  solutions[i] = res
end
simres2 = ConvergenceSimulation(solutions)

alg = "CrankNicholson"
solutions = cell(N)
for i = 1:N
  Δx = 1//2^(topΔx-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg)
  solutions[i] = res
end
simres3 = ConvergenceSimulation(solutions)

convplot_fullΔx(simres,titleStr="")
convplot_fullΔx(simres2,titleStr="")
convplot_fullΔx(simres3,titleStr="Dx Convergence Plots",savefile="dxconv.svg")
