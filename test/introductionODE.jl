# Introduction to the ODE Solvers

using DifferentialEquations

### Defining a problem
#First we define a problem type by giving it the equation and the initial condition
"""Example problem with solution ``u(t)=u₀*exp(α*t)``"""
function linearODEExample(;α=1,u₀=1/2)
  f(u,t) = α*u
  sol(u₀,t) = u₀*exp(α*t)
  return(ODEProblem(f,u₀,sol=sol))
end
prob = linearODEExample()
Δt = 1//2^(4) #The initial timestepping size
T = 1 # The final time

### Solve and plot
println("Solve and Plot")
sol =solve(prob::ODEProblem,Δt,T,fullSave=true,alg="Euler")
plot(sol,plottrue=true)
#Use Plots.jl's gui() command to display the plot.
Plots.gui()
#Shown is both the true solution and the approximated solution.

### Extras
#We can choose a better method as follows:
sol =solve(prob::ODEProblem,1//2^(4),1,fullSave=true,alg="ExplicitRK")
plot(sol,plottrue=true)
Plots.gui()
#This is pretty exact, so we can turn on adaptive timestepping to solve in less steps
sol =solve(prob::ODEProblem,1//2^(4),1,fullSave=true,alg="ExplicitRK",adaptive=true)
plot(sol,plottrue=true)
Plots.gui()
#More features can be designated via keyword arguments. Please see the manual for details.
true
