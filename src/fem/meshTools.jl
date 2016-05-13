"""
FEMmesh

"""
type FEMmesh
  node
  elem
  bdNode
  freeNode
  bdEdge
  isBdNode
  isBdElem
  bdFlag
  totalEdge
  area
  Dirichlet
  Neumann
  Robin
  N::Int
  NT::Int
  Δx
  Δt
  T::Number
  numIters::Int
  μ
  ν
  evolutionEq
  function FEMmesh(node,elem,Δx,Δt,T,bdType)
    N = size(node,1); NT = size(elem,1);
    totalEdge = [elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]]

    #Compute the area of each element
    ve = Array{Float64}(size(node[elem[:,3],:])...,3)
    ## Compute vedge, edge as a vector, and area of each element
    ve[:,:,1] = node[elem[:,3],:]-node[elem[:,2],:]
    ve[:,:,2] = node[elem[:,1],:]-node[elem[:,3],:]
    ve[:,:,3] = node[elem[:,2],:]-node[elem[:,1],:]
    area = 0.5*abs(-ve[:,1,3].*ve[:,2,2]+ve[:,2,3].*ve[:,1,2])

    #Boundary Conditions
    bdNode,bdEdge,isBdNode,isBdElem = findboundary(elem)
    bdFlag = setboundary(node::AbstractArray,elem::AbstractArray,bdType)
    Dirichlet = totalEdge[vec(bdFlag .== 1),:]
    Neumann = totalEdge[vec(bdFlag .== 2),:]
    Robin = totalEdge[vec(bdFlag .== 3),:]
    isBdNode = falses(N,1)
    isBdNode[Dirichlet] = true
    bdNode = find(isBdNode)
    freeNode = find(!isBdNode)
    if Δt != 0
      numIters = round(Int64,T/Δt)
    else
      numIters = 0
    end
    new(node,elem,bdNode,freeNode,bdEdge,isBdNode,isBdElem,bdFlag,totalEdge,area,Dirichlet,Neumann,Robin,N,NT,Δx,Δt,T,numIters,CFLμ(Δt,Δx),CFLν(Δt,Δx),T!=0)
  end
  FEMmesh(node,elem,Δx,bdType)=FEMmesh(node,elem,Δx,0,0,bdType)
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
meshgrid(vx)

Computes an (x,y)-grid from the vectors (vx,vx).
For more information, see the MATLAB documentation.
"""
meshgrid(v::AbstractVector) = meshgrid(v, v)

"""
meshgrid(vx,vy)

Computes an (x,y)-grid from the vectors (vx,vy).
For more information, see the MATLAB documentation.
"""
function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repmat(vx, m, 1), repmat(vy, 1, n))
end

"""
meshgrid(vx,vy,vz)

Computes an (x,y,z)-grid from the vectors (vx,vy,vz).
For more information, see the MATLAB documentation.
"""
function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T},
                     vz::AbstractVector{T})
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

"""
notime_squaremesh(square,Δx,bdType)

Computes the (node,elem) square mesh for the square
with the chosen Δx and boundary settings.

###Example
```juila
square=[0 1 0 1] #Unit Square
Δx=.25
notime_squaremesh(square,Δx,"Dirichlet")
```
"""
function notime_squaremesh(square,Δx,bdType)
  node,elem = fem_squaremesh(square,Δx)
  return(FEMmesh(node,elem,Δx,bdType))
end

"""
parabolic_squaremesh(square,Δx,Δt,T,bdType)

Computes the (node,elem) x [0,T] parabolic square mesh
for the square with the chosen Δx and boundary settings
and with the constant time intervals Δt.

###Example
```juila
square=[0 1 0 1] #Unit Square
Δx=.25; Δt=.25;T=2
parabolic_squaremesh(square,Δx,Δt,T,"Dirichlet")
```
"""
function parabolic_squaremesh(square,Δx,Δt,T,bdType)
  node,elem = fem_squaremesh(square,Δx)
  return(FEMmesh(node,elem,Δx,Δt,T,bdType))
end
