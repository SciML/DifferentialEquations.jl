##Finite Element Method Introduction

using DifferentialEquations, Plots

### Setup

#First we have to generate a mesh
Δx = 1//2^(3)
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet)
#Then we define our problem. To do this, you only need to define the equations for the PDE

f(x)=sin(2π.*x[:,1]).*cos(2π.*x[:,2])
prob = PoissonProblem(f)
#=
Here we have the true solution and the true gradient `Du`. The solvers will
automatically use these to calculate errors. Now we generate the problem type:
=#

### Solving

#To solve an FEMProblem, we only need to pass the solver the mesh and the problem
sol = solve(fem_mesh,prob)
#=
The solver picks the dispatch for the Poisson example, and by default solvers via
A direct solve with \. Returned is a solver object with all the knowledge of the solution.
=#

### Plotting

#=
The plotting abilities are given by Plots.jl. Since DifferentialEquations.jl
defines recipes for the solution objects, we can plot the default plot via:
=#
gr()
TEST_PLOT && plot(sol::FEMSolution) #To save the plot, use savefig("plot.png") or "plot.pdf", etc.

### Test Results

#The test will only pass if the calculated L2 error is below 1e-4.
true
