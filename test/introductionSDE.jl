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
Δt = 1//2^(4) #The initial timestepping size
T = 1 # The final time

### Solve and plot
#We can plot using the classic Euler-Maruyama algorithm as follows:
sol =solve(prob::SDEProblem,Δt,T,fullSave=true,alg="EM")
plot(sol,plottrue=true)
#Use Plots.jl's gui() command to display the plot.
gui()

### Extras
#We can choose a better method as follows:
sol =solve(prob::SDEProblem,Δt,T,fullSave=true,alg="SRI")
plot(sol,plottrue=true)
gui()
