
<a id='Internal-Finite-Element-Tools-1'></a>

# Internal Finite Element Tools


<a id='General-1'></a>

## General

<a id='DifferentialEquations' href='#DifferentialEquations'>#</a>
**`DifferentialEquations`** &mdash; *Module*.



###DifferentialEquations

This is a package for solving numerically solving differential equations in Julia by Chris Rackauckas. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include stochastic ordinary differential equations (SODEs or SDEs), stochastic partial differential equations (SPDEs), partial differential equations (with both finite difference and finite element methods), and differential delay equations. For ordinary differential equation solvers, see [ODE.jl](https://github.com/JuliaLang/ODE.jl)

This package is for efficient and parallel implementations of research-level algorithms, many of which are quite recent. These algorithms aim to be optimized for HPC applications, including the use of GPUs, Xeon Phis, and multi-node parallelism. With the easy to use plot/convergence testing algorithms, this package also provides a good sandbox for developing novel numerical schemes.


<a id='Mesh-Tools-1'></a>

## Mesh Tools

<a id='DifferentialEquations.meshgrid' href='#DifferentialEquations.meshgrid'>#</a>
**`DifferentialEquations.meshgrid`** &mdash; *Function*.



meshgrid(vx,vy,vz)

Computes an (x,y,z)-grid from the vectors (vx,vy,vz). For more information, see the MATLAB documentation.

meshgrid(vx,vy)

Computes an (x,y)-grid from the vectors (vx,vy). For more information, see the MATLAB documentation.

meshgrid(vx)

Computes an (x,y)-grid from the vectors (vx,vx). For more information, see the MATLAB documentation.

<a id='DifferentialEquations.CFLν' href='#DifferentialEquations.CFLν'>#</a>
**`DifferentialEquations.CFLν`** &mdash; *Function*.



CFLν(Δt,Δx)

Computes the CFL-condition ν= Δt/Δx

<a id='DifferentialEquations.CFLμ' href='#DifferentialEquations.CFLμ'>#</a>
**`DifferentialEquations.CFLμ`** &mdash; *Function*.



CFLμ(Δt,Δx)

Computes the CFL-condition μ= Δt/(Δx*Δx)


<a id='Solver-Tools-1'></a>

## Solver Tools

<a id='DifferentialEquations.gradbasis' href='#DifferentialEquations.gradbasis'>#</a>
**`DifferentialEquations.gradbasis`** &mdash; *Function*.



gradbasis(node,elem)

Returns the gradient of the barycentric basis elements.

<a id='DifferentialEquations.quadfbasis' href='#DifferentialEquations.quadfbasis'>#</a>
**`DifferentialEquations.quadfbasis`** &mdash; *Function*.



quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear;gNquadorder=2)

<a id='DifferentialEquations.quadpts' href='#DifferentialEquations.quadpts'>#</a>
**`DifferentialEquations.quadpts`** &mdash; *Function*.



quadpts(order)

Returns the quadrature points and weights for and order ### quadrature in 2D.

<a id='DifferentialEquations.accumarray' href='#DifferentialEquations.accumarray'>#</a>
**`DifferentialEquations.accumarray`** &mdash; *Function*.



accumarray(subs, val, sz=(maximum(subs),))

See MATLAB's documentation for more details.

<a id='DifferentialEquations.assemblematrix' href='#DifferentialEquations.assemblematrix'>#</a>
**`DifferentialEquations.assemblematrix`** &mdash; *Function*.



assemblematrix(node,elem;lumpflag=false,K=[])

Assembles the stiffness matrix A as an approximation to Δ on the finite element mesh (node,elem). Also generates the mass matrix M. If lumpflag=true, then the mass matrix is lumped resulting in a diagonal mass matrix. Specify a diffusion constant along the nodes via K.

### Returns

A = Stiffness Matrix M = Mass Matrix area = A vector of the calculated areas for each element.

assemblematrix(FEMmesh::FEMmesh;lumpflag=false,K=[])

Assembles the stiffness matrix A as an approximation to Δ on the finite element mesh (node,elem). Also generates the mass matrix M. If lumpflag=true, then the mass matrix is lumped resulting in a diagonal mass matrix. Specify a diffusion constant along the nodes via K.

### Returns

A = Stiffness Matrix M = Mass Matrix area = A vector of the calculated areas for each element.

<a id='DifferentialEquations.gradu' href='#DifferentialEquations.gradu'>#</a>
**`DifferentialEquations.gradu`** &mdash; *Function*.



gradu(node,elem,u,Dlambda=[])

Estimates the gradient of u on the mesh (node,elem)


<a id='Error-Tools-1'></a>

## Error Tools

<a id='DifferentialEquations.getH1error' href='#DifferentialEquations.getH1error'>#</a>
**`DifferentialEquations.getH1error`** &mdash; *Function*.



function getH1error(node,elem,Du,uh,K=[],quadOrder=[])

getH1error(femMesh::FEMmesh,Du,u)

Estimates the H1 error between uexact and uh on the mesh (node,elem). It reads the mesh to estimate the element type and uses this to choose a quadrature order unless specified. If K is specified then it is the diffusion coefficient matrix.

<a id='DifferentialEquations.getL2error' href='#DifferentialEquations.getL2error'>#</a>
**`DifferentialEquations.getL2error`** &mdash; *Function*.



getL2error(node,elem,uexact,uh,quadOrder=[])

getL2error(femMesh::FEMmesh,sol,u)

Estimates the L2 error between uexact and uh on the mesh (node,elem). It reads the mesh to estimate the element type and uses this to choose a quadrature order unless specified.

