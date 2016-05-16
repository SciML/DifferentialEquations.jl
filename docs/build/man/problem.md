
<a id='Defining-a-Problem-1'></a>

# Defining a Problem


Below are the definitions of the types which specify problems. Some general notes are:


  * (x,t) vs (x,y,t): Mathematically one normally specifies equations in 2D as `f(x,y,t)`. However, in this code we use `x` as a vector. Thus you can think of `x`=`x[:,1]` and `y`=`x[:,2]`. Thus input equations are of the form `f(x,t)` no matter the dimension. If time is not included in the problem (for example, a Poisson equation problem), then we use `f(x)`. An example is the equation `u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)` would be specified as `sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)`.
  * Linearity: If the equation has linear term, they are specified with functions `f(x,t)`. If it is nonlinear, it is specified with functions `f(u,x,t)`. The boundary conditions are always `(x,t)`
  * Stochastic: By default the equation is deterministic. For each equation, one can specify a σ term which adds a stochastic `σ(u,x,t)dW_t` term to the equation (or with `σ(x,t)dW_t` if linear, must match `f`). `dW_t` corresponds to the type of noise which is chosen. By default this is space-time Gaussian white noise.


<a id='Poisson-Equation-Problem-1'></a>

## Poisson Equation Problem

<a id='DifferentialEquations.PoissonProblem' href='#DifferentialEquations.PoissonProblem'>#</a>
**`DifferentialEquations.PoissonProblem`** &mdash; *Type*.



PoissonProblem

Wraps the data that define a 2D linear Poisson equation problem:

```math
-Δu = f
```

with bounday conditions `gD` on the Dirichlet boundary and gN on the Neumann boundary. Linearity is determined by whether the forcing function `f` is a function of two variables (x,t) or three (u,x,t) (with x=[:,1] and y=[:,2]).

If they keyword `σ` is given, then this wraps the data that define a 2D stochastic heat equation

```math
-Δu = f + σdW
```

###Constructors

PoissonProblem(f,sol,Du): Defines the Dirichlet problem with solution `sol`, solution gradient `Du = [u_x,u_y]`, and forcing function `f`

PoissonProblem(u₀,f): Defines the problem with initial value `u₀` (as a function) and f. If your initial data is a vector, wrap it as u₀(x) = vector.

Note: If all functions are of (x,t), then the program assumes it's linear. Write your functions using x = x[:,1] and y = x[:,2].  Use f=f(u,x,t) and σ=σ(u,x,t) (if specified) for nonlinear problems (with the boundary conditions still (x,t))

###Keyword Arguments

  * `gD` = Dirichlet boundary function

  * `gN` = Neumann boundary function

  * `σ` = The function which multiplies the noise `dW`. By default `σ` is 0.

  * `noiseType` = A string which specifies the type of noise to be generated. By default `noiseType` is "White" for Gaussian Spacetime White Noise.


<a id='Heat-Equation-Problem-1'></a>

## Heat Equation Problem

<a id='DifferentialEquations.HeatProblem' href='#DifferentialEquations.HeatProblem'>#</a>
**`DifferentialEquations.HeatProblem`** &mdash; *Type*.



HeatProblem

Wraps the data that define a 2D heat equation problem:

```math
u_t = Δu + f
```

with bounday conditions `gD` on the Dirichlet boundary and gN on the Neumann boundary. Linearity is determined by whether the forcing function `f` is a function of two variables (x,t) or three (u,x,t) (with x=[:,1] and y=[:,2]).

If they keyword `σ` is given, then this wraps the data that define a 2D stochastic heat equation

```math
u_t = Δu + f + σdW_t
```

###Constructors

  * `HeatProblem(sol,Du,f)`: Defines the Dirichlet problem with solution `sol`, solution gradient `Du = [u_x,u_y]`, and the forcing function `f`.

  * `HeatProblem(u₀,f)`: Defines the problem with initial value `u₀` (as a function) and `f`. If your initial data is a vector, wrap it as u₀(x) = vector.

Note: If all functions are of (x,t), then the program assumes it's linear. Write your functions using x = x[:,1] and y = x[:,2].  Use f=f(u,x,t) and σ=σ(u,x,t) (if specified) for nonlinear problems (with the boundary conditions still (x,t))

###Keyword Arguments

  * `gD` = Dirichlet boundary function

  * `gN` = Neumann boundary function

  * `σ` = The function which multiplies the noise dW. By default σ is 0.

  * `noiseType` = A string which specifies the type of noise to be generated. By default noiseType is "White" for Gaussian Spacetime White Noise.


<a id='Example-Problems-1'></a>

## Example Problems

<a id='DifferentialEquations.poissonProblemExample_wave' href='#DifferentialEquations.poissonProblemExample_wave'>#</a>
**`DifferentialEquations.poissonProblemExample_wave`** &mdash; *Function*.



Example problem with solution: `u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)`

<a id='DifferentialEquations.poissonProblemExample_noisyWave' href='#DifferentialEquations.poissonProblemExample_noisyWave'>#</a>
**`DifferentialEquations.poissonProblemExample_noisyWave`** &mdash; *Function*.



Example problem with deterministic solution: `u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)`

<a id='DifferentialEquations.poissonProblemExample_birthdeath' href='#DifferentialEquations.poissonProblemExample_birthdeath'>#</a>
**`DifferentialEquations.poissonProblemExample_birthdeath`** &mdash; *Function*.



Example problem for nonlinear Poisson equation. Uses `f(u)=1-u/2`.

<a id='DifferentialEquations.heatProblemExample_diffuse' href='#DifferentialEquations.heatProblemExample_diffuse'>#</a>
**`DifferentialEquations.heatProblemExample_diffuse`** &mdash; *Function*.



Example problem with solution: `u(x,y,t)=exp(-10((x-.5).^2 + (y-.5).^2 )-t)`

<a id='DifferentialEquations.heatProblemExample_pure' href='#DifferentialEquations.heatProblemExample_pure'>#</a>
**`DifferentialEquations.heatProblemExample_pure`** &mdash; *Function*.



Example problem which starts with 1 at (0.5,0.5) and solves with `f=gD=0`

<a id='DifferentialEquations.heatProblemExample_moving' href='#DifferentialEquations.heatProblemExample_moving'>#</a>
**`DifferentialEquations.heatProblemExample_moving`** &mdash; *Function*.



Example problem with solution: `u(x,y,t)=0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x-t+0.5).^2 + (y-t+0.5).^2))`

<a id='DifferentialEquations.heatProblemExample_birthdeath' href='#DifferentialEquations.heatProblemExample_birthdeath'>#</a>
**`DifferentialEquations.heatProblemExample_birthdeath`** &mdash; *Function*.



Example problem which starts with 0 and solves with `f(u)=1-u/2`

<a id='DifferentialEquations.heatProblemExample_stochasticbirthdeath' href='#DifferentialEquations.heatProblemExample_stochasticbirthdeath'>#</a>
**`DifferentialEquations.heatProblemExample_stochasticbirthdeath`** &mdash; *Function*.



Example problem which starts with 0 and solves with `f(u)=1-u/2` with noise `σ(u)=10u^2`


<a id='Related-Functions-1'></a>

## Related Functions

<a id='DifferentialEquations.PdeProblem' href='#DifferentialEquations.PdeProblem'>#</a>
**`DifferentialEquations.PdeProblem`** &mdash; *Type*.



PdeProblem: Defines PDE problems via its internal functions

