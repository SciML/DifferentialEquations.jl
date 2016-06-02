# Ordinary Differential Equation (ODE) Example

In this example we will solve the equation

```math
uₜ = f(u,t)dt
```

where ``f(u,t)=αu``. We know via Calculus that the solution to this equation is
``u(t)=u₀*exp(α*t)``. To solve this numerically, we define a problem type by
giving it the equation and the initial condition:

```julia
"""Example problem with solution ``u(t)=u₀*exp(α*t)``"""
function linearODEExample(;α=1,u₀=1/2)
  f(u,t) = α*u
  sol(u₀,t) = u₀*exp(α*t)
  return(ODEProblem(f,u₀,sol=sol))
end
```

Then we setup some parameters:

```julia
Δt = 1//2^(4) #The initial timestepping size
T = 1 # The final time
```

We then send these items to the solver.

```julia
sol =solve(prob::ODEProblem,Δt,T,fullSave=true,alg="Euler")
```

Plotting commands are provided via a recipe to Plots.jl. To plot the solution
object, simply call plot:

```julia
plot(sol,plottrue=true)
#Use Plots.jl's gui() command to display the plot.
Plots.gui()
#Shown is both the true solution and the approximated solution.
```

More keyword arguments can be found in the Plots.jl documentation. When we plot
this solution, we see that it is off from the true solution. We can choose a
better algorithm by specifying:

```julia
sol =solve(prob::ODEProblem,1//2^(4),1,fullSave=true,alg="ExplicitRK")
plot(sol,plottrue=true)
Plots.gui()
```

The `"ExplicitRK"` algorithms are general Runge-Kutta solvers. It defaults to
Dormand-Prince 4/5, the same solver as MATLAB's `ode45`. Please see the solver
documentation for more algorithms.

Notice that this solution tracks the true solution really well. Thus we can
solve the problem in less timesteps by turning on adaptive timestepping. To
do so, you simply pass a keyword argument:

```julia
sol =solve(prob::ODEProblem,1//2^(4),1,fullSave=true,alg="ExplicitRK",adaptive=true)
plot(sol,plottrue=true)
Plots.gui()
```
