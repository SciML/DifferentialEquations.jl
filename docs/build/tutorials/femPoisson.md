
<a id='Poisson-Equation-Finite-Element-Method-Example-1'></a>

# Poisson Equation Finite Element Method Example


In this example we will solve the Poisson Equation Δu=f. The code for this example can be found in [tests/introductionExample.jl](tests/introductionExample.jl). For our example, we will take the linear equation where `f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])`. For this equation we know that solution is `u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)` with gradient `Du(x,y) = [cos(2*pi.*x).*cos(2*pi.*y)./(4*pi) -sin(2π.*x).*sin(2π.*y)./(4π)]`. Thus, we define a PoissonProblem as follows:


```julia
"Example problem with solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  gN(x) = 0
  isLinear = true
  return(PoissonProblem(f,sol,Du,gN,isLinear))
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

