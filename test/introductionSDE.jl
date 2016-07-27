# Introduction to the SDE Solvers

using DifferentialEquations, Plots
import DifferentialEquations.linearSDEExample #Ignore this

srand(100)
### Defining a problem
#First we define a problem type by giving it the equation and the initial condition
"""Example problem with solution ``u(t,W)=u₀*exp((α-(β^2)/2)*t+β*W)``"""
function linearSDEExample(;α=1,β=1,u₀=1/2)
  f(u,t) = α*u
  σ(u,t) = β*u
  sol(u₀,t,W) = u₀*exp((α-(β^2)/2)*t+β*W)
  return(SDEProblem(f,σ,u₀,sol=sol))
end
prob = linearSDEExample()
Δt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.

### Solve and plot
#We can solve using the classic Euler-Maruyama algorithm:
sol =solve(prob::SDEProblem,tspan,Δt=Δt,fullSave=true,alg=:EM)
plot(sol,plottrue=true)
#Use Plots.jl's gui() command to display the plot.
gui()

### Extras
#We can choose a better method as follows:
sol =solve(prob::SDEProblem,tspan,Δt=Δt,fullSave=true,alg=:SRI)
plot(sol,plottrue=true)
gui()
