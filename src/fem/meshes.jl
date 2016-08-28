"""
FEMmesh

Holds the information describing a finite element mesh. For information on how (node,elem)
can be interpreted as a mesh describing a geometry, see the mesh specification documentation.

### Fields

* `node`: The nodes in the (node,elem) structure.
* `elem`: The elements in the (node,elem) structure.
* `bdnode`: Vector of indices for the boundary nodes.
* `freenode`: Vector of indices for the free (non-dirichlet bound) nodes.
* `bdedge`: Indices of the edges in totaledge which are on the boundary.
* `is_bdnode`: Boolean which is true for nodes on the boundary.
* `is_bdelem`: Boolean which is true for elements on the boundary.
* `bdflag`: Flag which describes the type of boundary condition. 1=> dirichlet,
2=>neumann, 3=>robin.
* `totaledge`: Vector of the edges.
* `area`: Vector which is the area for each element.
* `dirichlet`: Indices for the nodes on the boundary which have a dirichlet boundary condition.
* `neumann`: Indices for the nodes on the boundary which have a neumann boundary condition.
* `robin`: Indices for the nodes on the boundary which have a robin boundary condition.
* `N::Int`: The number of nodes.
* `NT`::Int: The number of triangles (elements).
* `Δx`: The spatial discretization size. If non-uniform, this is the average.
* `Δt`: The time discretization size. If adaptive, this is the initial.
* `T`::Number: The end time.
* `numiters`::Int: The number of iterations to go from 0 to T using Δt.
* `μ`: The CFL μ stability parameter.
* `ν`: The CFL ν stability parameter.
* `evolutionEq`: True for a mesh which has non-trivial time components.

"""
type FEMmesh <: Mesh
  node
  elem
  bdnode
  freenode
  bdedge
  is_bdnode
  is_bdelem
  bdflag
  totaledge
  area
  dirichlet
  neumann
  robin
  N::Int
  NT::Int
  Δx
  Δt
  T::Number
  numiters::Int
  μ
  ν
  evolutionEq
  function FEMmesh(node,elem,Δx,Δt,T,bdtype)
    N = size(node,1); NT = size(elem,1);
    totaledge = [elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]]

    #Compute the area of each element
    ve = Array{eltype(node)}(size(node[elem[:,3],:])...,3)
    ## Compute vedge, edge as a vector, and area of each element
    ve[:,:,1] = node[elem[:,3],:]-node[elem[:,2],:]
    ve[:,:,2] = node[elem[:,1],:]-node[elem[:,3],:]
    ve[:,:,3] = node[elem[:,2],:]-node[elem[:,1],:]
    area = 0.5*abs(-ve[:,1,3].*ve[:,2,2]+ve[:,2,3].*ve[:,1,2])

    #Boundary Conditions
    bdnode,bdedge,is_bdnode,is_bdelem = findboundary(elem)
    bdflag = setboundary(node::AbstractArray,elem::AbstractArray,bdtype)
    dirichlet = totaledge[vec(bdflag .== 1),:]
    neumann = totaledge[vec(bdflag .== 2),:]
    robin = totaledge[vec(bdflag .== 3),:]
    is_bdnode = falses(N,1)
    is_bdnode[dirichlet] = true
    bdnode = find(is_bdnode)
    freenode = find(!is_bdnode)
    if Δt != 0
      numiters = round(Int64,T/Δt)
    else
      numiters = 0
    end
    new(node,elem,bdnode,freenode,bdedge,is_bdnode,is_bdelem,bdflag,totaledge,area,dirichlet,neumann,robin,N,NT,Δx,Δt,T,numiters,CFLμ(Δt,Δx),CFLν(Δt,Δx),T!=0)
  end
  FEMmesh(node,elem,Δx,bdtype)=FEMmesh(node,elem,Δx,0,0,bdtype)
end

"""
SimpleMesh

Holds the information describing a finite element mesh. For information on how (node,elem)
can be interpreted as a mesh describing a geometry, see [Programming of Finite
Element Methods by Long Chen](http://www.math.uci.edu/~chenlong/226/Ch3FEMCode.pdf).

### Fields

* `node`: The nodes in the (node,elem) structure.
* `elem`: The elements in the (node,elem) structure.
"""
type SimpleMesh <: Mesh
  node
  elem
  function SimpleMesh(node,elem)
    return(new(node,elem))
  end
end


"""
CFLμ(Δt,Δx)

Computes the CFL-condition μ= Δt/(Δx*Δx)
"""
CFLμ(Δt,Δx)=Δt/(Δx*Δx)

"""
CFLν(Δt,Δx)

Computes the CFL-condition ν= Δt/Δx
"""
CFLν(Δt,Δx)=Δt/Δx

"""
fem_squaremesh(square,h)

Returns the grid in the iFEM form of the two arrays (node,elem)
"""
function fem_squaremesh(square,h)
  x0 = square[1]; x1= square[2];
  y0 = square[3]; y1= square[4];
  x,y = meshgrid(x0:h:x1,y0:h:y1)
  node = [x[:] y[:]];

  ni = size(x,1); # number of rows
  N = size(node,1);
  t2nidxMap = 1:N-ni;
  topNode = ni:ni:N-ni;
  t2nidxMap = deleteat!(collect(t2nidxMap),collect(topNode));
  k = t2nidxMap;
  elem = [k+ni k+ni+1 k ; k+1 k k+ni+1];
  return(node,elem)
end

"""
notime_squaremesh(square,Δx,bdtype)

Computes the (node,elem) square mesh for the square
with the chosen Δx and boundary settings.

###Example

```julia
square=[0 1 0 1] #Unit Square
Δx=.25
notime_squaremesh(square,Δx,"dirichlet")
```
"""
function notime_squaremesh(square,Δx,bdtype)
  node,elem = fem_squaremesh(square,Δx)
  return(FEMmesh(node,elem,Δx,bdtype))
end

"""
parabolic_squaremesh(square,Δx,Δt,T,bdtype)

Computes the (node,elem) x [0,T] parabolic square mesh
for the square with the chosen Δx and boundary settings
and with the constant time intervals Δt.

###Example
```julia
square=[0 1 0 1] #Unit Square
Δx=.25; Δt=.25;T=2
parabolic_squaremesh(square,Δx,Δt,T,:dirichlet)
```
"""
function parabolic_squaremesh(square,Δx,Δt,T,bdtype)
  node,elem = fem_squaremesh(square,Δx)
  return(FEMmesh(node,elem,Δx,Δt,T,bdtype))
end

type FDMMesh
  Δxs::AbstractArray
  mins::AbstractArray
  maxs::AbstractArray
  grids::Tuple
  dims::Int
  gridSize::Tuple
  square::Bool
  function FDMMesh(Δxs::AbstractArray;mins=[0;0],maxs=[1;1],buildMesh=true)
    if length(mins)!=length(maxs) || length(Δxs)!=length(mins) error("DimensionMismatch") end
    dims = length(mins)
    if buildMesh && dims == 2
      grids = meshgrid(mins[1]:Δxs[1]:maxs[1],mins[2]:Δxs[2]:maxs[2])
      gridSize = size(grids[1])
    elseif buildMesh && dims == 3
      grids = meshgrid(mins[1]:Δxs[1]:maxs[1],mins[2]:Δxs[2]:maxs[2],mins[3]:Δxs[3]:maxs[3])
      gridSize = size(grids[1])
    else
      grids = nothing
      gridSize = (maxs-mins)./Δxs
    end
    square = minimum(map((x)->(x == gridSize[1]),gridSize))
    new(Δxs,mins,maxs,grids,dims,gridSize,square)
  end
  FDMMesh(Δx::Number;mins=[0;0],maxs=[1;1],buildMesh=true) = FDMMesh(Δx*ones(mins),mins=mins,maxs=maxs,buildMesh=buildMesh)
end

"""
size(mesh::FDMMesh)

Returns gridSize.
"""
size(mesh::FDMMesh) = mesh.gridSize
