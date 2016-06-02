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

### Other Algorithms

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

### Systems of Equations

We can also solve systems of equations. DifferentialEquations.jl can handle any
size problem, so instead of showing it for a vector, let's let u be a matrix!
To do this, we simply need to have u₀ be a matrix, and define f such that it
takes in a matrix and outputs a matrix. We can define a matrix of linear ODEs
as follows:

```julia
"""Example problem of 8 linear ODEs (as a 4x2 matrix) with
solution ``u(t)=exp(α.*t)`` and random initial conditions"""
function twoDimlinearODEExample(;α=ones(4,2),u₀=rand(4,2).*ones(4,2)/2)
  f(u,t) = α.*u
  sol(u₀,t) = u₀.*exp(α.*t)
  return(ODEProblem(f,u₀,sol=sol))
end
prob = twoDimlinearODEExample()
```

Here our ODE is on a 4x2 matrix. Since we are using .\*, this is 8 independent
ODEs, but you can do whatever you want. To solve the ODE, we do the same steps
as before.

```julia
sol =solve(prob::ODEProblem,1//2^(4),1,fullSave=true,alg="ExplicitRK")
plot(sol,plottrue=true)
Plots.gui()
```

Notice now we have 8 solutions and 8 true solutions, but since we used the high
order method, the true solutions are covered by the approximations.
