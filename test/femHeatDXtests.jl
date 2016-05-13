######
##FEM Heat Δx Convergence Tests
######
using DifferentialEquations

T = 1
#Travis CI Test Setting
#Not good plots, but quick for unit tests
Δt = 1//2^(6) #Small Δt for Euler stability, but takes long
N = 2
topΔx = 3

if !isdefined(:testState) #Don't plot during test
  # Convergence Test Configuration
  # Use this setup to get good plots
  Δt = 1//2^(14) #Small Δt for Euler stability, but takes long
  N = 4
  topΔx = 7
end
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

if !isdefined(:testState) #Don't plot during test
  convplot_fullΔx(simres,titleStr="")
  convplot_fullΔx(simres2,titleStr="")
  convplot_fullΔx(simres3,titleStr="Dx Convergence Plots",savefile="dxconv.svg")
end
