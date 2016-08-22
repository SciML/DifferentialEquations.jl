# Stochastic Differential Equation (SDE) Example

This tutorial will introduce you to the functionality for solving SDE. Other
introductions can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/examples).

In this example we will solve the equation

```math
du = f(u,t)dt + Σσᵢ(u,t)dWⁱ
```

where ``f(u,t)=αu`` and ``σ(u,t)=βu``. We know via Stochastic Calculus that the
solution to this equation is ``u(t,W)=u₀\exp((α-\frac{β^2}{2})t+βW)``. To solve this
numerically, we define a problem type by giving it the equation and the initial
condition:

```julia
using DifferentialEquations
α=1
β=1
u₀=1/2
f(u,t) = α*u
σ(u,t) = β*u
Δt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.
```

For reference, let's also give the `SDEProblem` the analytical solution. Note that
each of the problem types allow for this, but it's always optional. This can be
a good way to judge how accurate the algorithms are, or is used to test convergence
of the algorithms for methods developers. Thus we define the problem object with:

```julia
analytic(u₀,t,W) = u₀*exp((α-(β^2)/2)*t+β*W)
prob = SDEProblem(f,σ,u₀,analytic=analytic)
```

and then we pass this information to the solver and plot:

```julia
#We can plot using the classic Euler-Maruyama algorithm as follows:
sol =solve(prob::SDEProblem,tspan,Δt=Δt,alg=:EM)
plot(sol,plot_analytic=true)
#Use Plots.jl's gui() command to display the plot.
Plots.gui()
```

<img src="https://raw.githubusercontent.com/ChrisRackauckas/DifferentialEquations.jl/master/examples/plots/introSDEplot.png" width="750" align="middle"  />

We can choose a higher-order solver for a more accurate result:

```julia
#We can choose a better method as follows:
sol =solve(prob::SDEProblem,tspan,Δt=Δt,alg=:SRIW1Optimized)
plot(sol,plot_analytic=true)
Plots.gui()
```

<img src="https://raw.githubusercontent.com/ChrisRackauckas/DifferentialEquations.jl/master/examples/plots/introSDEplotSRI.png" width="750" align="middle"  />
