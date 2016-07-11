"""
findboundary(elem,bdFlag=[])

findboundary(femMesh::FEMmesh,bdFlag=[])

Finds elements which are on the boundary of the domain. If bdFlag is given,
then those indices are added as nodes for a Dirichlet boundary condition (useful
for creating cracks and other cutouts of domains).

### Returns
bdNode = Vector of indices for bdNode. Using node[:,bdNode] returns boundary nodes.

bdEdge = Vector of indices for boundary edges.

isBdNode = Vector of booleans size N which donotes which are on the boundary

isBdElem = Vector of booleans size NT which denotes which are on the boundary

"""
function findboundary(elem::AbstractArray;bdFlag=[])
  N = round(Int,maximum(elem))
  elem = round(Int,elem)
  nv = size(elem,2)
  if nv == 3 # triangle
      totalEdge = [elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]]
  elseif nv == 4
      totalEdge = [elem[:,[1,2]]; elem[:,[2,3]]; elem[:,[3,4]]; elem[:,[4,1]]]
  end
  if !isempty(bdFlag)
      Dirichlet = totalEdge[(vec(bdFlag) .== 1),:]
      isBdNode = falses(N)
      isBdNode[vec(Dirichlet)] = true
      bdNode = find(isBdNode)
      bdEdge = totalEdge[(vec(bdFlag) .== 2) | (vec(bdFlag) .== 3),:]
  else
      totalEdge = sort(totalEdge,2)
      edgeMatrix = sparse(totalEdge[:,1],totalEdge[:,2],1)
      i,j = ind2sub(size(edgeMatrix),find(x->x==1,edgeMatrix))
      bdEdge = [i';j']'
      isBdNode = falses(N)
      isBdNode[bdEdge] = true
      bdNode = find(isBdNode)
  end
  isBdElem = isBdNode[elem[:,1]] | isBdNode[elem[:,2]] | isBdNode[elem[:,3]]
  return(bdNode,bdEdge,isBdNode,isBdElem)
end

findboundary(femMesh::Mesh,bdFlag=[]) = findboundary(femMesh.elem,bdFlag=bdFlag)

"""
setboundary(node::AbstractArray,elem::AbstractArray,bdType)

setboundary(femMesh::FEMmesh,bdType)

Takes in the femMesh and creates an array bdFlag which denotes the boundary types.
1 stands for Dirichlet, 2 for Neumann, 3 for Robin.
"""
function setboundary(node::AbstractArray,elem::AbstractArray,bdType)
  ## Find boundary edges
  nv = size(elem,2)
  if nv == 3 # triangles
      totalEdge = sort([elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]],2)
  elseif nv == 4 # quadrilateral
      totalEdge = sort([elem[:,[1,2]]; elem[:,[2,3]]; elem[:,[3,4]]; elem[:,[4,1]]],2)
  end
  ne = nv # number of edges in one element
  Neall = size(totalEdge,1)
  NT = size(elem,1)
  edge = unique(totalEdge,1)
  totalEdge = sort(totalEdge,2)
  edgeMatrix = sparse(round(Int,totalEdge[:,1]),round(Int,totalEdge[:,2]),1)
  i,j = ind2sub(size(edgeMatrix),find(x->x==1,edgeMatrix))
  bdEdge = [i';j']'
  bdEdgeidx = zeros(Int64,size(bdEdge,1))
  for i = 1:size(bdEdge,1)
    bdEdgeidx[i] = find(all(totalEdge .== bdEdge[i,:]', 2))[1] #Find the edge in totalEdge and save index
  end

  bdFlag = zeros(Int8,Neall)
  ## Set up boundary edges
  #nVarargin = size(varargin,2)
  #if (nVarargin==1)
  bdType = findbdtype(bdType)
  bdFlag[bdEdgeidx] = bdType
  #end
  #=
  if (nVarargin>=2)
      for i=1:nVarargin/2
          bdType = findbdtype(varargin{2*i-1})
          expr = varargin{2*i}
          if strcmp(expr,"all")
              bdFlag(bdEdgeidx) = bdType
          else
             x = (node(allEdge(bdEdgeidx,1),1) + node(allEdge(bdEdgeidx,2),1))/2
             y = (node(allEdge(bdEdgeidx,1),2) + node(allEdge(bdEdgeidx,2),2))/2
             idx = eval(expr)
             bdFlag(bdEdgeidx(idx)) = bdType
          end
      end
  end
  =#
  bdFlag = reshape(bdFlag,NT,ne)
  return(bdFlag)
end

setboundary(femMesh::Mesh,bdType) = setboundary(femMesh.node,femMesh.elem,bdType)

function findbdtype(bdstr)
        if bdstr=="Dirichlet"
            bdType = 1
        elseif bdstr=="Neumann"
            bdType = 2
        elseif bdstr=="Robin"
            bdType = 3
        #=
        elseif bdstr=="ABC" # absorbing boundary condition for wave-type equations
            bdType = 4
        =#
    end
    return(bdType)
end
