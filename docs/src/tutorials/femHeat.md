# Heat Equation Finite Element Method Example

In this example we will solve the heat equation ``u_t=Δu+f``. To do this, we define
a HeatProblem which contains the function `f` and the boundary conditions. We
specify one as follows:

```julia
"Example problem which starts with 0 and solves with f(u)=1-.1u"
function heatProblemExample_birthdeath()
  gD(x,t) = zeros(size(x,1))
  f(u,x,t)  = ones(size(x,1)) - .5u
  u0(x) = zeros(size(x,1))
  gN(x,t) = 0
  isLinear = false
  return(HeatProblem(u0,f,gD,gN,isLinear))
end
pdeProb = heatProblemExample_birthdeath()
```

Here the equation we chose was nonlinear since `f` depends on the variable `u`.
Thus we specify f=f(u,x,t) and set isLinear = false. If `f` did not depend on
u, then we would specify f=f(x,t) and `isLinear = true`. `gD` specifies the condition
on the Dirichlet part of the boundary and `gN` specifies the condition on the
Neumann part of the boundary. `u0` specifies the initial condition. These together
give a HeatProblem object which holds everything which specifies a Heat Equation Problem.

We then generate a mesh. We will solve the equation on the parabolic cylinder
[0,1]^2 x [0,1]. You can think of this as the cube, or at every time point from 0
to 1, the domain is the unit square. To generate this mesh, we use the command

```julia
T = 1
Δx = 1//2^(3)
Δt = 1//2^(7)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
```  

We then call the appropriate solver

```julia
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
```

Here we have chosen to use the Euler algorithm to solve the equation. Other algorithms
and their descriptions can be found in the solvers part of the documentation.
