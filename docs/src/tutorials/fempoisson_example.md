# Poisson Equation Finite Element Method Example

This tutorial will introduce you to the functionality for solving a PDE. Other
introductions can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/examples).

In this example we will solve the Poisson Equation ``Δu=f``. For our example, we will take the linear equation where ``f(x,y) = \sin(2πx)\cos(2πy)``. For this equation we know that solution is ``u(x,y,t)= \sin(2πx)\cos(2πy)/(8π^2)`` with gradient ``Du(x,y) = [\cos(2πx)\cos(2πy)/(4π) -\sin(2πx)\sin(2πy)/(4π)]``. Thus, we define a PoissonProblem as follows:

```julia
f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
gD(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
prob = PoissonProblem(f,gD)
```

Here we chose the dirichlet boundary condition `gD` to give the theoretical solution.  Other example problems can be found in [src/examples/exampleProblems.jl](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/src/premades/premade_problems.jl). To solve this problem, we first have to generate a mesh. Here we will simply generate a mesh of triangles on the square [0,1]x[0,1] with Δx=2^(-5). To do so, we use the code:

```julia
Δx = 1//2^(5)
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet)
```

Note that by specifying :dirichlet our boundary conditions is set on all boundaries to dirichlet. This gives an FEMmesh object which stores a finite element mesh in the same layout as [iFEM](http://www.math.uci.edu/~chenlong/programming.html). Notice this code shows that the package supports the use of rationals in meshes. Other numbers such as floating point and integers can be used as well. Finally, to solve the equation we use

```julia
sol = solve(fem_mesh,pdeProb)
```

solve takes in a mesh and a PoissonProblem and uses the solver to compute the solution. Here the solver was chosen to be GMRES. Other solvers can be found in the documentation. This returns a FEMSolution object which holds data about the solution, such as the solution values (u). To plot the solution, we use the command

```julia
plot(sol::FEMSolution)
Plots.gui()
```

Here is the plot shown against the analytical solution to show the accuracy:

<img src="https://raw.githubusercontent.com/ChrisRackauckas/DifferentialEquations.jl/master/examples/plots/introductionExample.png" width="750" align="middle"  />
