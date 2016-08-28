# Meshes

## Mesh Specification

Finite element meshes are specified in the (node,elem) structure due to Long Chen.
For the standard elements used in this package, we describe a geometric figure by
a triangulation. The nodes are the vertices of the triangle and the elements are the
triangles themselves. These are encoded as follows:

* Row ``i`` of node is an ``(x,y)`` (or ``(x,y,z)``) pair which specifies the coordinates
  of the ``i``th node.
* Row ``j`` of elem are the indices of the nodes which make the triangle. Thus in
  2D each row has three numbers.

For example, to know the ``(x,y)`` locations of the vertices of triangle ``j``, we
would see that `node[elem[j,i],:]` are the ``(x,y)`` locations of the ``i``th vertex
for ``i=1,2,3``.

For more information, please see [Programming of Finite
Element Methods by Long Chen](http://www.math.uci.edu/~chenlong/226/Ch3FEMCode.pdf).

## Mesh Type

```@docs
DifferentialEquations.FEMmesh
DifferentialEquations.SimpleMesh
DifferentialEquations.Mesh
```

## Mesh Generation Functions

```@docs
DifferentialEquations.findboundary
DifferentialEquations.setboundary
DifferentialEquations.fem_squaremesh
DifferentialEquations.notime_squaremesh
DifferentialEquations.parabolic_squaremesh
Base.size(::DifferentialEquations.FDMMesh)
```

## Example Meshes

```@docs
DifferentialEquations.meshExample_bunny
DifferentialEquations.meshExample_flowpastcylindermesh
DifferentialEquations.meshExample_lakemesh
DifferentialEquations.meshExample_Lshapemesh
DifferentialEquations.meshExample_Lshapeunstructure
DifferentialEquations.meshExample_oilpump
DifferentialEquations.meshExample_wavymesh
DifferentialEquations.meshExample_wavyperturbmesh
```

## Plot Functions

The plot functionality is provided by a Plots.jl recipe. What is plotted is a
"trisurf" of the mesh. To plot a mesh, simply use:

```julia
plot(mesh::Mesh)
```

All of the functionality (keyword arguments) provided by Plots.jl are able to
be used in this command. Please see the Plots.jl documentation for more information.
