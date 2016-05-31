######
##FEM Poisson Nonlinear Tests
######
using DifferentialEquations

Δx = 1//2^(3)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Neumann")
prob = poissonProblemExample_birthdeath()

sol = solve(femMesh::FEMmesh,prob::PoissonProblem)

plot(sol,plottrue=false,zlim=(0,2))

#Returns true if computed solution is homogenous near 2
maximum(abs(sol.u - 2))< 1e-9
