
<a id='Mesh-Generation-1'></a>

# Mesh Generation


<a id='Mesh-Specification-1'></a>

## Mesh Specification


<a id='Mesh-Type-1'></a>

## Mesh Type

<a id='DifferentialEquations.FEMmesh' href='#DifferentialEquations.FEMmesh'>#</a>
**`DifferentialEquations.FEMmesh`** &mdash; *Type*.



FEMmesh


<a id='Mesh-Generation-Functions-1'></a>

## Mesh Generation Functions

<a id='DifferentialEquations.findboundary' href='#DifferentialEquations.findboundary'>#</a>
**`DifferentialEquations.findboundary`** &mdash; *Function*.



findboundary(elem,bdFlag=[])

findboundary(femMesh::FEMmesh,bdFlag=[])

Finds elements which are on the boundary of the domain. If bdFlag is given, then those indices are added as nodes for a Dirichlet boundary condition (useful for creating cracks and other cutouts of domains).

### Returns

bdNode = Vector of indices for bdNode. Using node[:,bdNode] returns boundary nodes.

bdEdge = Vector of indices for boundary edges.

isBdNode = Vector of booleans size N which donotes which are on the boundary

isBdElem = Vector of booleans size NT which denotes which are on the boundary

<a id='DifferentialEquations.setboundary' href='#DifferentialEquations.setboundary'>#</a>
**`DifferentialEquations.setboundary`** &mdash; *Function*.



setboundary(node::AbstractArray,elem::AbstractArray,bdType)

setboundary(femMesh::FEMmesh,bdType)

Takes in the femMesh and creates an array bdFlag which denotes the boundary types. 1 stands for Dirichlet, 2 for Neumann, 3 for Robin. 

<a id='DifferentialEquations.fem_squaremesh' href='#DifferentialEquations.fem_squaremesh'>#</a>
**`DifferentialEquations.fem_squaremesh`** &mdash; *Function*.



fem_squaremesh(square,h)

Returns the grid in the iFEM form of the two arrays (node,elem)

<a id='DifferentialEquations.notime_squaremesh' href='#DifferentialEquations.notime_squaremesh'>#</a>
**`DifferentialEquations.notime_squaremesh`** &mdash; *Function*.



notime_squaremesh(square,Δx,bdType)

Computes the (node,elem) square mesh for the square with the chosen Δx and boundary settings.

###Example `square=[0 1 0 1] #Unit Square` `Δx=.25` `notime_squaremesh(square,Δx,"Dirichlet")`

<a id='DifferentialEquations.parabolic_squaremesh' href='#DifferentialEquations.parabolic_squaremesh'>#</a>
**`DifferentialEquations.parabolic_squaremesh`** &mdash; *Function*.



parabolic_squaremesh(square,Δx,Δt,T,bdType)

Computes the (node,elem) x [0,T] parabolic square mesh for the square with the chosen Δx and boundary settings and with the constant time intervals Δt.

###Example `square=[0 1 0 1] #Unit Square` `Δx=.25; Δt=.25;T=2` `parabolic_squaremesh(square,Δx,Δt,T,"Dirichlet")`

