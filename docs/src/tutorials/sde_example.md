# Stochastic Differential Equation (SDE) Example

This tutorial will introduce you to the functionality for solving SDE. Other
introductions can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/examples).

In this example we will solve the equation

```math
du = f(t,u)dt + Σσᵢ(t,u)dWⁱ
```

where ``f(t,u)=αu`` and ``σ(t,u)=βu``. We know via Stochastic Calculus that the
solution to this equation is ``u(t,W)=u₀\exp((α-\frac{β^2}{2})t+βW)``. To solve this
numerically, we define a problem type by giving it the equation and the initial
condition:

```julia
using DifferentialEquations
α=1
β=1
u₀=1/2
f(t,u) = α*u
σ(t,u) = β*u
Δt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.
```

For reference, let's also give the `SDEProblem` the analytical solution. Note that
each of the problem types allow for this, but it's always optional. This can be
a good way to judge how accurate the algorithms are, or is used to test convergence
of the algorithms for methods developers. Thus we define the problem object with:

```julia
analytic(t,u₀,W) = u₀*exp((α-(β^2)/2)*t+β*W)
prob = SDEProblem(f,σ,u₀,analytic=analytic)
```

and then we pass this information to the solver and plot:

```julia
#We can plot using the classic Euler-Maruyama algorithm as follows:
sol =solve(prob::SDEProblem,tspan,Δt=Δt,alg=:EM)
plot(sol,plot_analytic=true)
```

![SDE Solution](https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/examples/plots/introSDEplot.png)

We can choose a higher-order solver for a more accurate result:

```julia
sol =solve(prob::SDEProblem,tspan,Δt=Δt,alg=:SRIW1Optimized)
plot(sol,plot_analytic=true)
```

![Better SDE Solution](https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/examples/plots/introSDEplotSRI.png)
