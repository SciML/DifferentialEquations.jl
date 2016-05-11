######
##FEM Poisson Δx Convergence Tests
######
using DiffEq

N = 4
topΔx = 7
pdeProb = poissonProblemExample_wave()

solutions = cell(N)
for i = 1:N
  Δx = 1//2^(topΔx-i)
  femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
  res = fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem)
  solutions[i] = res
end

simres = ConvergenceSimulation(solutions)
convplot_fullΔx(simres,titleStr="Poisson Δx Convergence")
