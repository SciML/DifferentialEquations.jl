######
##FEM Poisson Nonlinear System Tests
######
using DifferentialEquations

Δx = 1//2^(1)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Neumann")
prob = poissonProblemExample_birthdeathsystem()

sol = solve(femMesh::FEMmesh,prob::PoissonProblem)

plot(sol,plottrue=false,zlim=(0,2))

#Returns true if computed solution is homogenous near 2
bool1 = maximum(abs(sol.u .- [2 1]))< 1e-8

### Harder system

prob = poissonProblemExample_birthdeathinteractingsystem()

sol = solve(femMesh::FEMmesh,prob::PoissonProblem)

plot(sol,plottrue=false,zlim=(0,2),cbar=false)

bool2 = maximum(abs(sol.u .- [2 1]))< 1e-8

bool1 && bool2
