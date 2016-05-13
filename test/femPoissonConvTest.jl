######
##FEM Poisson Δx Convergence Tests
######
using DifferentialEquations,LaTeXStrings

N = 2 # 2 for testing, use 4 for good graph
topΔx = 4 # 4 for testing, use 7 for good graph
pdeProb = poissonProblemExample_wave()

solutions = cell(N)
for i = 1:N #Loop through Δx's, solve Poisson, and save to solutions
  Δx = 1//2^(topΔx-i)
  femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
  res = fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem)
  solutions[i] = res
end

#Construct ConvergenceSimulation Object
simres = ConvergenceSimulation(solutions)

#Plot Result
dxstring = L"\Delta x"
if !isdefined(:testState) #Don't plot during test
  convplot_fullΔx(simres,titleStr="Poisson $dxstring Convergence",savefile="plot.png")
end
