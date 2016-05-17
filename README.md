# DifferentialEquations.jl

[![Join the chat at https://gitter.im/ChrisRackauckas/DifferentialEquations.jl](https://badges.gitter.im/ChrisRackauckas/DifferentialEquations.jl.svg)](https://gitter.im/ChrisRackauckas/DifferentialEquations.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![Build Status](https://travis-ci.org/ChrisRackauckas/DifferentialEquations.jl.svg?branch=master)](https://travis-ci.org/ChrisRackauckas/DifferentialEquations.jl) [![Build status](https://ci.appveyor.com/api/projects/status/032otj4kh462tq2l/branch/master?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/differentialequations-jl/branch/master) [![Coverage Status](https://coveralls.io/repos/github/ChrisRackauckas/DifferentialEquations.jl/badge.svg?branch=master)](https://coveralls.io/github/ChrisRackauckas/DifferentialEquations.jl?branch=master) [![codecov](https://codecov.io/gh/ChrisRackauckas/DifferentialEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ChrisRackauckas/DifferentialEquations.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ChrisRackauckas.github.io/DifferentialEquations.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://ChrisRackauckas.github.io/DifferentialEquations.jl/latest)

This is a package for solving numerically solving differential equations in Julia by Chris Rackauckas. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include stochastic ordinary differential equations (SODEs or SDEs), stochastic partial differential equations (SPDEs), partial differential equations (with both finite difference and finite element methods), and differential delay equations. For ordinary differential equation solvers, see [ODE.jl](https://github.com/JuliaLang/ODE.jl)

This package is for efficient and parallel implementations of research-level algorithms, many of which are quite recent. These algorithms aim to be optimized for HPC applications, including the use of GPUs, Xeon Phis, and multi-node parallelism. With the easy to use plot/convergence testing algorithms, this package also provides a good sandbox for developing novel numerical schemes. Since this package is designed for long computations, one of the features of this package is the existence of tools for inspecting a long calculation. These include optional printing and, if the user is using Juno, a progress meter (with time estimates once implemented on Juno's end).

Currently, finite element solvers for the (Stochastic) Poisson and Heat Equations are supplied. These functions take in a problem specification (with an option for adding stochasticity) and generate solutions to the PDEs. Mesh generation tools currently only work for squares, though the solvers will work on general meshes if they are provided. The mesh layout follows the format of [iFEM](http://www.math.uci.edu/~chenlong/programming.html) and many subroutines in the finite element solver are based off of iFEM algorithms.

If you have any questions, or just want to chat about solvers/using the package, please feel free to message me in the Gitter channel. For bug reports, feature requests, etc., please submit an issue.

# Using the package

To install the package, use the following command inside the Julia REPL:
```julia
Pkg.add("DifferentialEquations")
```

For all of the latest features, switch to the master branch via:

```julia
Pkg.checkout("DifferentialEquations")
```

To load the package, use the command

```julia
using DifferentialEquations
```

To understand the package in more detail, check out the examples codes in [test/](test/).
Note that for many of the examples in the test folder, you may wish to run them at
lower Δx or Δt. These values were taken to be large in order make unit tests run faster!
For the most up to date information, please contact me [via the repository Gitter](https://gitter.im/ChrisRackauckas/DifferentialEquations.jl)
or [read the latest documentation](http://chrisrackauckas.github.io/DifferentialEquations.jl/latest/)

## Poisson Equation Finite Element Method Example

In this example we will solve the Poisson Equation Δu=f. The code for this example can be found in [test/introductionExample.jl](test/introductionExample.jl). For our example, we will take the linear equation where `f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])`. For this equation we know that solution is `u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)` with gradient `Du(x,y) = [cos(2*pi.*x).*cos(2*pi.*y)./(4*pi) -sin(2π.*x).*sin(2π.*y)./(4π)]`. Thus, we define a PoissonProblem as follows:

```julia
"Example problem with solution: u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  return(PoissonProblem(f,sol,Du))
end
pdeProb = poissonProblemExample_wave()
```

Note that in this case since the solution is known, the Dirichlet boundary condition `gD` is automatically set to match the true solution. The code for other example problems can be found in [src/examples/exampleProblems.jl](src/examples/exampleProblems.jl). To solve this problem, we first have to generate a mesh. Here we will simply generate a mesh of triangles on the square [0,1]x[0,1] with Δx=2^(-5). To do so, we use the code:

```julia
Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
```

Note that by specifying "Dirichlet" our boundary conditions is set on all boundaries to Dirichlet. This gives an FEMmesh object which stores a finite element mesh in the same layout as [iFEM](http://www.math.uci.edu/~chenlong/programming.html). Notice this code shows that the package supports the use of rationals in meshes. Other numbers such as floating point and integers can be used as well. Finally, to solve the equation we use

```julia
res = fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem,solver="GMRES")
```

fem_solvepoisson takes in a mesh and a PoissonProblem and uses the solver to compute the solution. Here the solver was chosen to be GMRES. Other solvers can be found in the documentation. This reurns a FEMSolution object which holds data about the solution, such as the solution values (u), the true solution (uTrue), error estimates, etc. To plot the solution, we use the command

```julia
solplot(res,savefile="introductionExample.png")
```

This gives us the following plot:

<img src="/src/examples/introductionExample.png" width="750" align="middle"  />

### Finite Element Stochastic Poisson Equation

We can solve the same PDE as a stochastic PDE, -Δu=f+gdW, with additive space-time white noise by specifying the problem as:

```julia
"Example problem with deterministic solution: u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_noisyWave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  σ(x) = 5 #Additive noise, a big amount!
  return(PoissonProblem(f,sol,Du,σ=σ))
end
```

and using the same solving commands as shown in [femStochasticPoissonSolo.jl](/src/femStochasticPoissonSolo.jl). This gives the following plot:

<img src="/src/examples/introductionStochasticExample.png" width="750" align="middle" />

### Finite Element Stochastic Heat Equation

The last equation we will solve in this introductory example will be a nonlinear stochastic heat equation u_t=Δu+f+gdW with forcing function `f(u)=.5-u`, noise function `g(u)=100u^2` and
initial condition `u0=0`. We would expect this system to rise towards the deterministic steady state `u=2` (but stay in mean a bit below it due to 1st order "Milstein" effects), gaining more noise as it increases. This is specified as follows:

```julia
"Example problem which starts with 0 and solves with f(u)=1-.1u"
function heatProblemExample_stochasticbirthdeath()
  f(u,x,t)  = ones(size(x,1)) - .5u
  u₀(x) = zeros(size(x,1))
  σ(u,x,t) = 100u.^2
  return(HeatProblem(u₀,f,σ=σ,stochastic=stochastic))
end
```

As shown in [femStochasticHeatAnimationTest.jl](/src/femStochasticHeatAnimationTest.jl), we use the following code create an animation of the solution:

```julia
T = 5
Δx = 1//2^(4)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
pdeProb = heatProblemExample_stochasticbirthdeath()

res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler",fullSave=true)
solplot_animation(res::FEMSolution;zlim=(0,2),cbar=false)
```

<img src="/src/examples/stochasticHeatAnimation.gif" width="750" align="middle" />

# Supported Equations
* (Stochastic) PDE Solvers
  * Finite Element Solvers
    * Linear Poisson Equation
    * Semi-linear Poisson Equation
    * Linear Heat Equation
    * Semi-linear Heat Equation (aka Reaction-Diffusion Equation)

# Roadmap
* SODE Solvers
  * Euler-Maruyama
  * Milstein
  * Rossler-SRK
  * Adaptive-SRK
* (Stochastic) PDE Solvers
  * Finite difference solvers:
    * Semi-linear Heat Equation (Reaction-Diffusion Equation)
    * Semi-linear Poisson Equation
    * Wave Equation
    * Transport Equation
    * Stokes Equation
    * Implicit Integration Factor (IIF) Maruyama
    * Implicit Integration Factor (IIF) Milstein
