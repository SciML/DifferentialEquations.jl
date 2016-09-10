# Stochastic Finite Element Examples

This tutorial will introduce you to the functionality for solving SPDEs. Other
introductions can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/examples).

For most PDE problem types, we can additionally specify them as a stochastic
problem by giving the appropriate optional arguments to the constructor. These
arguments are a function σ which is the function multiplied to the Brownian
increments ``dW``, and stochastic, a boolean which we put as true for when the equation
is stochastic. Another keyword that is optional is noisetype which specifies the
type of noise (the "color" of the noise). By default this is Gaussian (Space-time)
White Noise.

The following examples show how to change the tutorial problems into stochastic problems.

## Finite Element Stochastic Poisson Equation

We can solve the same PDE as in the Poisson Tutorial except as the stochastic PDE,
 ``-Δu=f+gdW``, with additive space-time white noise by specifying the problem as:

```julia
f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
σ(x) = 5 #Additive noise
prob = PoissonProblem(f,σ=σ)
solve(prob)
```

This gives the following plot (with adding the deterministic solution from the previous example):

![Stochastic Poisson Example Solution](https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/examples/plots/introductionStochasticExample.png)

## Finite Element Stochastic Heat Equation

This will solve a nonlinear stochastic heat equation u_t=Δu+f+gdW with forcing function `f(u)=.5-u`,
noise function `g(u)=100u^2` and initial condition `u0=0`. We would expect this system
to rise towards the deterministic steady state `u=2` (but stay in mean a bit below
it due to 1st order "Milstein" effects), gaining more noise as it increases.
This is specified as follows:

```julia
f(t,x,u)  = ones(size(x,1)) - .5u
u₀(x) = zeros(size(x,1))
σ(t,x,u) = 1u.^2
prob = HeatProblem(u₀,f,σ=σ)
```

We use the following code create an animation of the solution:

```julia
T = 5
Δx = 1//2^(3)
Δt = 1//2^(11)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)

sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler,save_timeseries=true,solver=:LU)
animate(sol::FEMSolution;zlim=(0,3),cbar=false)
```

![Stochastic Heat Solution](https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/examples/plots/stochasticHeatAnimation.gif)
