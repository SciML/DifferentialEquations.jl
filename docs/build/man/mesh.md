
<a id='Meshes-1'></a>

# Meshes


<a id='Mesh-Specification-1'></a>

## Mesh Specification


Finite element meshes are specified in the (node,elem) structure due to Long Chen. For the standard elements used in this package, we describe a geometric figure by a triangulation. The nodes are the vertices of the triangle and the elements are the triangles themselves. These are encoded as follows:


  * Row `i` of node is an `(x,y)` (or `(x,y,z)`) pair which specifies the coordinates of the `i`th node.
  * Row `j` of elem are the indices of the nodes which make the triangle. Thus in 2D each row has three numbers.


For example, to know the `(x,y)` locations of the vertices of triangle `j`, we would see that `node[elem[j,i],:]` are the `(x,y)` locations of the `i`th vertex for `i=1,2,3`.


For more information, please see [Programming of Finite Element Methods by Long Chen](http://www.math.uci.edu/~chenlong/226/Ch3FEMCode.pdf).


<a id='Mesh-Type-1'></a>

## Mesh Type

<a id='DifferentialEquations.FEMmesh' href='#DifferentialEquations.FEMmesh'>#</a>
**`DifferentialEquations.FEMmesh`** &mdash; *Type*.



FEMmesh

Holds the information describing a finite element mesh. For information on how (node,elem) can be interpreted as a mesh describing a geometry, see the mesh specification documentation.

**Fields**

  * `node`: The nodes in the (node,elem) structure.
  * `elem`: The elements in the (node,elem) structure.
  * `bdNode`: Vector of indices for the boundary nodes.
  * `freeNode`: Vector of indices for the free (non-Dirichlet bound) nodes.
  * `bdEdge`: Indices of the edges in totalEdge which are on the boundary.
  * `isBdNode`: Boolean which is true for nodes on the boundary.
  * `isBdElem`: Boolean which is true for elements on the boundary.
  * `bdFlag`: Flag which describes the type of boundary condition. 1=> Dirichlet, 2=>Neumann, 3=>Robin.
  * `totalEdge`: Vector of the edges.
  * `area`: Vector which is the area for each element.
  * `Dirichlet`: Indices for the nodes on the boundary which have a Dirichlet boundary condition.
  * `Neumann`: Indices for the nodes on the boundary which have a Neumann boundary condition.
  * `Robin`: Indices for the nodes on the boundary which have a Robin boundary condition.
  * `N::Int`: The number of nodes.
  * `NT`::Int: The number of triangles (elements).
  * `Δx`: The spatial discretization size. If non-uniform, this is the average.
  * `Δt`: The time discretization size. If adaptive, this is the initial.
  * `T`::Number: The end time.
  * `numIters`::Int: The number of iterations to go from 0 to T using Δt.
  * `μ`: The CFL μ stability parameter.
  * `ν`: The CFL ν stability parameter.
  * `evolutionEq`: True for a mesh which has non-trivial time components.

<a id='DifferentialEquations.SimpleMesh' href='#DifferentialEquations.SimpleMesh'>#</a>
**`DifferentialEquations.SimpleMesh`** &mdash; *Type*.



SimpleMesh

Holds the information describing a finite element mesh. For information on how (node,elem) can be interpreted as a mesh describing a geometry, see [Programming of Finite Element Methods by Long Chen](http://www.math.uci.edu/~chenlong/226/Ch3FEMCode.pdf).

**Fields**

  * `node`: The nodes in the (node,elem) structure.
  * `elem`: The elements in the (node,elem) structure.

<a id='DifferentialEquations.Mesh' href='#DifferentialEquations.Mesh'>#</a>
**`DifferentialEquations.Mesh`** &mdash; *Type*.



Mesh: An abstract type which holds a (node,elem) pair and other information for a mesh


<a id='Mesh-Generation-Functions-1'></a>

## Mesh Generation Functions

<a id='DifferentialEquations.findboundary' href='#DifferentialEquations.findboundary'>#</a>
**`DifferentialEquations.findboundary`** &mdash; *Function*.



findboundary(elem,bdFlag=[])

findboundary(femMesh::FEMmesh,bdFlag=[])

Finds elements which are on the boundary of the domain. If bdFlag is given, then those indices are added as nodes for a Dirichlet boundary condition (useful for creating cracks and other cutouts of domains).

**Returns**

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

###Example

```julia
square=[0 1 0 1] #Unit Square
Δx=.25
notime_squaremesh(square,Δx,"Dirichlet")
```

<a id='DifferentialEquations.parabolic_squaremesh' href='#DifferentialEquations.parabolic_squaremesh'>#</a>
**`DifferentialEquations.parabolic_squaremesh`** &mdash; *Function*.



parabolic_squaremesh(square,Δx,Δt,T,bdType)

Computes the (node,elem) x [0,T] parabolic square mesh for the square with the chosen Δx and boundary settings and with the constant time intervals Δt.

###Example

```julia
square=[0 1 0 1] #Unit Square
Δx=.25; Δt=.25;T=2
parabolic_squaremesh(square,Δx,Δt,T,"Dirichlet")
```


<a id='Example-Meshes-1'></a>

## Example Meshes

<a id='DifferentialEquations.meshExample_bunny' href='#DifferentialEquations.meshExample_bunny'>#</a>
**`DifferentialEquations.meshExample_bunny`** &mdash; *Function*.



meshExample_bunny() : Returns a 3D SimpleMesh.

<a id='DifferentialEquations.meshExample_flowpastcylindermesh' href='#DifferentialEquations.meshExample_flowpastcylindermesh'>#</a>
**`DifferentialEquations.meshExample_flowpastcylindermesh`** &mdash; *Function*.



meshExample_flowpastcylindermesh() : Returns a 2D SimpleMesh.

<a id='DifferentialEquations.meshExample_lakemesh' href='#DifferentialEquations.meshExample_lakemesh'>#</a>
**`DifferentialEquations.meshExample_lakemesh`** &mdash; *Function*.



meshExample_lakemesh() : Returns a 2D SimpleMesh.

<a id='DifferentialEquations.meshExample_Lshapemesh' href='#DifferentialEquations.meshExample_Lshapemesh'>#</a>
**`DifferentialEquations.meshExample_Lshapemesh`** &mdash; *Function*.



meshExample_Lshapemesh() : Returns a 2D SimpleMesh.

<a id='DifferentialEquations.meshExample_Lshapeunstructure' href='#DifferentialEquations.meshExample_Lshapeunstructure'>#</a>
**`DifferentialEquations.meshExample_Lshapeunstructure`** &mdash; *Function*.



meshExample_Lshapeunstructure() : Returns a 2D SimpleMesh.

<a id='DifferentialEquations.meshExample_oilpump' href='#DifferentialEquations.meshExample_oilpump'>#</a>
**`DifferentialEquations.meshExample_oilpump`** &mdash; *Function*.



meshExample_oilpump() : Returns a 3D SimpleMesh.

<a id='DifferentialEquations.meshExample_wavymesh' href='#DifferentialEquations.meshExample_wavymesh'>#</a>
**`DifferentialEquations.meshExample_wavymesh`** &mdash; *Function*.



meshExample_wavymesh() : Returns a 2D SimpleMesh.

<a id='DifferentialEquations.meshExample_wavyperturbmesh' href='#DifferentialEquations.meshExample_wavyperturbmesh'>#</a>
**`DifferentialEquations.meshExample_wavyperturbmesh`** &mdash; *Function*.



meshExample_wavyperturbmesh() : Returns a 3D SimpleMesh.

