
<a id='Defining-a-Problem-1'></a>

# Defining a Problem


(x,t) vs (x,y,t) isLinear stochastic


<a id='Poisson-Equation-Problem-1'></a>

## Poisson Equation Problem

<a id='DifferentialEquations.PoissonProblem' href='#DifferentialEquations.PoissonProblem'>#</a>
**`DifferentialEquations.PoissonProblem`** &mdash; *Type*.



PoissonProblem

Wraps the data that define a 2D linear Poisson equation problem:

Δu = f(x,t)

#Constructors

PoissonProblem(f,sol,Du,gN,isLinear): Defines the Dirichlet problem with solution sol, solution gradient Du = [u_x,u_y], f, and Neumann boundary data gN,

PoissonProblem(u0,f,gD,gN,isLinear): Defines the problem with initial value u0 (as a function or vector), f, Dirichlet boundary function gD, and Neumann boundary function gN.

Note: If isLinear is true, then all functions must only be functions of (x). If isLinear is false, then f=f(u,x) and σ=σ(u,x) (if specified), while the other functions are only functions of (x).

#Keyword Arguments

σ = The function which multiplies the noise dW. By default σ is 0.

stochastic = A boolean which specifies if the problem is stochastic. By default stochastic is false.

noiseType = A string which specifies the type of noise to be generated. By default noiseType is "White" for Gaussian Spacetime White Noise.


<a id='Heat-Equation-Problem-1'></a>

## Heat Equation Problem

<a id='DifferentialEquations.HeatProblem' href='#DifferentialEquations.HeatProblem'>#</a>
**`DifferentialEquations.HeatProblem`** &mdash; *Type*.



`HeatProblem`

Wraps the data that define a 2D linear heat equation problem:

`u_t = Δu + f(x,t)`

#Constructors

HeatProblem(sol,Du,f,isLinear): Defines the Dirichlet problem with solution sol, solution gradient Du = [u_x,u_y], f, and a boolean which states whether the problem is linear (i.e. linear if f does not depend on u).

HeatProblem(u0,f,gD,gN,isLinear): Defines the problem with initial value u0 (as a function or vector), f, Dirichlet boundary function gD,  Neumann boundary function gN, and a boolean which states whether the problem is linear (i.e. linear if f does not depend on u).

Note: If isLinear is true, then all functions must only be functions of (x,t). If isLinear is false, then f=f(u,x,t) and σ=σ(u,x,t) (if specified), while the other functions are only functions of (x,t).

#Keyword Arguments

The constructors take the following keyword arguments:

σ = The function which multiplies the noise dW. By default σ is 0.

stochastic = A boolean which specifies if the problem is stochastic. By default stochastic is false.

noiseType = A string which specifies the type of noise to be generated. By default noiseType is "White" for Gaussian Spacetime White Noise.


<a id='Example-Problems-1'></a>

## Example Problems

<a id='DifferentialEquations.poissonProblemExample_wave' href='#DifferentialEquations.poissonProblemExample_wave'>#</a>
**`DifferentialEquations.poissonProblemExample_wave`** &mdash; *Function*.



Example problem with solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)

<a id='DifferentialEquations.poissonProblemExample_noisyWave' href='#DifferentialEquations.poissonProblemExample_noisyWave'>#</a>
**`DifferentialEquations.poissonProblemExample_noisyWave`** &mdash; *Function*.



Example problem with deterministic solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)

<a id='DifferentialEquations.heatProblemExample_diffuse' href='#DifferentialEquations.heatProblemExample_diffuse'>#</a>
**`DifferentialEquations.heatProblemExample_diffuse`** &mdash; *Function*.



Example problem with solution: u(x,y,t)=exp(-10((x-.5).^2 + (y-.5).^2 )-t)

<a id='DifferentialEquations.heatProblemExample_pure' href='#DifferentialEquations.heatProblemExample_pure'>#</a>
**`DifferentialEquations.heatProblemExample_pure`** &mdash; *Function*.



Example problem which starts with 1 at (0.5,0.5) and solves with f=gD=0

<a id='DifferentialEquations.heatProblemExample_moving' href='#DifferentialEquations.heatProblemExample_moving'>#</a>
**`DifferentialEquations.heatProblemExample_moving`** &mdash; *Function*.



Example problem with solution: u(x,y,t)=0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x-t+0.5).^2 + (y-t+0.5).^2))

<a id='DifferentialEquations.heatProblemExample_birthdeath' href='#DifferentialEquations.heatProblemExample_birthdeath'>#</a>
**`DifferentialEquations.heatProblemExample_birthdeath`** &mdash; *Function*.



Example problem which starts with 0 and solves with f(u)=1-.1u

<a id='DifferentialEquations.heatProblemExample_stochasticbirthdeath' href='#DifferentialEquations.heatProblemExample_stochasticbirthdeath'>#</a>
**`DifferentialEquations.heatProblemExample_stochasticbirthdeath`** &mdash; *Function*.



Example problem which starts with 0 and solves with f(u)=1-.1u


<a id='Related-Functions-1'></a>

## Related Functions

<a id='DifferentialEquations.PdeProblem' href='#DifferentialEquations.PdeProblem'>#</a>
**`DifferentialEquations.PdeProblem`** &mdash; *Type*.



PdeProblem: Defines PDE problems via its internal functions

