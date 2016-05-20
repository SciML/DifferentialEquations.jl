using DifferentialEquations

Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
pdeProb = poissonProblemExample_wave()

res = solve(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="CG")

mesh = SimpleMesh(femMesh.node,femMesh.elem)

true
