# Introduction to the ODE Solvers

using DifferentialEquations
import DifferentialEquations: linearODEExample, twoDimlinearODEExample # Ignore

### Defining a problem
#First we define a problem type by giving it the equation and the initial condition
"""Example problem with solution ``u(t)=u₀*exp(α*t)``"""
function linearODEExample(;α=1,u₀=1/2)
  f(u,t) = α*u
  sol(u₀,t) = u₀*exp(α*t)
  return(ODEProblem(f,u₀,sol=sol))
end
prob = linearODEExample()
Δt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.

#=
Note here we provided the true solution because it's known. However, the
true solution is optional and the solvers will work without it!
=#

### Solve and plot
println("Solve and Plot")
sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:Euler)
plot(sol,plottrue=true)
#Use Plots.jl's gui() command to display the plot.
Plots.gui()
#Shown is both the true solution and the approximated solution.

### Extras
#We can choose a better method as follows:
sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:ExplicitRK)
plot(sol,plottrue=true)
Plots.gui()
#This is pretty exact, so we can turn on adaptive timestepping to solve in less steps
sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:ExplicitRK,adaptive=true)
plot(sol,plottrue=true)
Plots.gui()
#More features can be designated via keyword arguments. Please see the manual for details.

## Multidimensional ODE
#=
Now we will solve a multidimensional ODE. DifferentialEquations.jl can handle any
size problem, so instead of showing it for a vector, let's let u be a matrix!
To do this, we simply need to have u₀ be a matrix, and define f such that it
takes in a matrix and outputs a matrix. We can define a matrix of linear ODEs
as follows:
=#
"""Example problem of 8 linear ODEs (as a 4x2 matrix) with solution ``u(t)=exp(α.*t)`` and random initial conditions"""
function twoDimlinearODEExample(;α=ones(4,2),u₀=rand(4,2).*ones(4,2)/2)
  f(u,t) = α.*u
  sol(u₀,t) = u₀.*exp(α.*t)
  return(ODEProblem(f,u₀,sol=sol))
end
prob = twoDimlinearODEExample()
#=
Here our ODE is on a 4x2 matrix. Since we are using .*, this is 8 independent
ODEs, but you can do whatever you want. To solve the ODE, we do the same steps
as before.
=#
sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:ExplicitRK)
plot(sol,plottrue=true)
Plots.gui()
#=
Notice now we have 8 solutions and 8 true solutions, but since we used the high
order method, the true solutions are covered by the approximations.
=#
true
