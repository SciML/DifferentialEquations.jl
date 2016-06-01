#Finite Element Method Introduction

using DifferentialEquations

### Setup

#First we have to generate a mesh
Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
#Then we define our problem. To do this, you only need to define the equations for the PDE
"Example problem with solution: ``u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)``"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  return(PoissonProblem(f,sol,Du))
end
#Here we have the true solution and the true gradient `Du`. The solvers will
#automatically use these to calculate errors. Now we generate the problem type:
pdeProb = poissonProblemExample_wave()

### Solving

#To solve an FEMProblem, we only need to pass the solver the mesh and the problem
sol = solve(femMesh,pdeProb)#,solver="CG") TODO fix CG
#The solver picks the dispatch for the Poisson example, and by default solvers via
#A direct solve with \. Returned is a solver object with all the knowledge of the solution.

### Plotting

#The plotting abilities are given by Plots.jl. Since DifferentialEquations.jl
#defines recipes for the solution objects, we can plot the default plot via:
plot(sol::FEMSolution,plottrue=false) #To save the plot, use savefig("plot.png") or "plot.pdf", etc.
#Note that if we set plottrue=true, we will plot the true solution alongside the
#approximated solution. Other arguments can be found via the Plots.jl documentation.

### Test Results

#The test will only pass if the calculated L2 error is below 1e-4.
sol.errors["L2"] < 1e-4
