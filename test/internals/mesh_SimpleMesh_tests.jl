using DifferentialEquations

Δx = 1//2^(5)
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet)
pdeProb = prob_poisson_wave

res = solve(fem_mesh::FEMmesh,pdeProb::PoissonProblem)#,solver=:CG) TODO Fix CG

mesh = SimpleMesh(fem_mesh.node,fem_mesh.elem)

true
