# Ordinary Differential Equation (ODE) Example

This tutorial will introduce you to the functionality for solving ODEs. Other
introductions can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/examples).

In this example we will solve the equation

```math
\frac{du}{dt} = f(t,u)
```

where ``f(t,u)=αu``. We know via Calculus that the solution to this equation is
``u(t)=u₀\exp(αt)``. To solve this numerically, we define a problem type by
giving it the equation and the initial condition:

```julia
using DifferentialEquations
α=1
u₀=1/2
f(t,u) = u
prob = ODEProblem(f,u₀)
```

Then we setup some parameters:

```julia
Δt = 1/2^(4) #The initial step size. It will automatically determined if not given.
tspan = [0,1] # The timespan. This is the default if not given.
```

We then send these items to the solver.

```julia
sol =solve(prob::ODEProblem,tspan,Δt=Δt,alg=:Euler)
```

To see what's in the solution object, we can print it:
```julia
print(sol)
#DifferentialEquations.ODESolution with 17 timesteps. No analytical solution is known.
#u: 1.3189642486832998
#t: [0.0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1.0]
#timeseries: [0.5,0.53125,0.564453,0.599731,0.637215,0.677041,0.719356,0.764315,0.812085,0.86284,0.916768,0.974066,1.03494,1.09963,1.16836,1.24138,1.31896]
```

We can access the 5th value of the solution with

```julia
sol[5]
#.637
```

or get the time of the 8th timestep by

```julia
sol.t[8]
#.438
```

The object that is returns by default acts as a continuous solution via an interpolation.
We can access the interpolated values by treating `sol` as a function, for example:

```julia
sol(0.45) # The value of the solution at t=0.45
```

For details on more finely controlling the output, see [the output specification manual page](http://juliadiffeq.github.io/DifferentialEquations.jl/latest/man/output_specification.html)

Plotting commands are provided via a recipe to Plots.jl. To plot the solution
object, simply call plot:

```julia
plot(sol)
#Use Plots.jl's gui() command to display the plot.
Plots.gui()
```

The plot function can be formatted using [the attributes available in Plots.jl](https://juliaplots.github.io/).
For more of an introduction to plotting solutions, [see the IJulia notebook](https://github.com/JuliaDiffEq/DifferentialEquations.jl/blob/master/examples/Formatting%20the%20Plots.ipynb).

### Other Algorithms

We can choose a better algorithm by specifying:

```julia
sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:ExplicitRK,adaptive=false)
plot(sol)
Plots.gui()
```

![Better ODE Solution](https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/examples/plots/introODEplot.png)


The `"ExplicitRK"` algorithms are general Runge-Kutta solvers. It defaults to
Dormand-Prince 4/5, the same solver as MATLAB's `ode45`. Please see the solver
documentation for more algorithms.

We can solve the problem in less timesteps by turning on adaptive timestepping. To
do so, you simply pass a keyword argument (note: this is true by default):

```julia
sol =solve(prob::ODEProblem,tspan,Δt=Δt,alg=:ExplicitRK,adaptive=true)
plot(sol)
Plots.gui()
```

![Adaptive ODE Solution](https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/examples/plots/adaptiveODEplot.png)

### Systems of Equations

We can also solve systems of equations. DifferentialEquations.jl can handle any
size problem, so instead of showing it for a vector, let's let u be a matrix!
To do this, we simply need to have u₀ be a matrix, and define f such that it
takes in a matrix and outputs a matrix. We can define a matrix of linear ODEs
as follows:

```julia
u₀=rand(4,2).*ones(4,2)/2
α=ones(4,2)
f(t,u) = α.*u
prob = ODEProblem(f,u₀)
```

Here our ODE is on a 4x2 matrix. Since we are using .\*, this is 8 independent
ODEs, but you can do whatever you want. To solve the ODE, we do the same steps
as before.

```julia
sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:ExplicitRK)
plot(sol)
Plots.gui()
```

![ODE System Solution](https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/examples/plots/multiODEplot.png)


### Defining Systems of Equations Eloquently Using @ode_define

To simplify your life, DifferentialEquations.jl provides the `@ode_define` macro
for "defining your ODE in pseudocode" and getting a function which is efficient
and runnable. For our example we will use [the Lorenz system](https://en.wikipedia.org/wiki/Lorenz_system).
The standard way to write this out in most mathematical programs is the following:

```julia
f = (t,u,du) -> begin
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
```

Here for more efficiency we plugged in the parameters. However, this does not
look like the pretty ``\LaTeX`` system we see on Wikipedia, and this might make it
harder to double-check that you defined the system correctly. Using the
`@ode_define` macro is much nicer:

```julia
g = @ode_define begin
  dx = σ*(y-x)
  dy = x*(ρ-z) - y
  dz = x*y - β*z
end σ=>10. ρ=>28. β=>(8/3)
```

DifferentialEquations.jl will automatically translate this to be exactly the
same as `f`. The result is more legible code with no performance loss.
