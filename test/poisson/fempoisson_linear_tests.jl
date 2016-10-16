##Finite Element Method Introduction

using DifferentialEquations, Plots

### Setup
Δx = 1//2^(3)
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet)

f = (x) -> sin.(2π.*x[:,1]).*cos.(2π.*x[:,2])
prob = PoissonProblem(f)

sol = solve(fem_mesh,prob)

sol = solve(fem_mesh,prob,solver=:CG)

sol = solve(fem_mesh,prob,solver=:GMRES)

TEST_PLOT && plot(sol::FEMSolution)

### Test Results

#The test will only pass if the calculated L2 error is below 1e-4.
true
