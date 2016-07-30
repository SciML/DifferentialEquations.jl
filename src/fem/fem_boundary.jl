"""
findboundary(elem,bdflag=[])

findboundary(fem_mesh::FEMmesh,bdflag=[])

Finds elements which are on the boundary of the domain. If bdflag is given,
then those indices are added as nodes for a dirichlet boundary condition (useful
for creating cracks and other cutouts of domains).

### Returns
bdnode = Vector of indices for bdnode. Using node[:,bdnode] returns boundary nodes.

bdedge = Vector of indices for boundary edges.

is_bdnode = Vector of booleans size N which donotes which are on the boundary

is_bdelem = Vector of booleans size NT which denotes which are on the boundary

"""
function findboundary(elem::AbstractArray;bdflag=[])
  N = round(Int,maximum(elem))
  elem = round(Int,elem)
  nv = size(elem,2)
  if nv == 3 # triangle
      totaledge = [elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]]
  elseif nv == 4
      totaledge = [elem[:,[1,2]]; elem[:,[2,3]]; elem[:,[3,4]]; elem[:,[4,1]]]
  end
  if !isempty(bdflag)
      dirichlet = totaledge[(vec(bdflag) .== 1),:]
      is_bdnode = falses(N)
      is_bdnode[vec(dirichlet)] = true
      bdnode = find(is_bdnode)
      bdedge = totaledge[(vec(bdflag) .== 2) | (vec(bdflag) .== 3),:]
  else
      totaledge = sort(totaledge,2)
      edge_matrix = sparse(totaledge[:,1],totaledge[:,2],1)
      i,j = ind2sub(size(edge_matrix),find(x->x==1,edge_matrix))
      bdedge = [i';j']'
      is_bdnode = falses(N)
      is_bdnode[bdedge] = true
      bdnode = find(is_bdnode)
  end
  is_bdelem = is_bdnode[elem[:,1]] | is_bdnode[elem[:,2]] | is_bdnode[elem[:,3]]
  return(bdnode,bdedge,is_bdnode,is_bdelem)
end

findboundary(fem_mesh::Mesh,bdflag=[]) = findboundary(fem_mesh.elem,bdflag=bdflag)

"""
setboundary(node::AbstractArray,elem::AbstractArray,bdtype)

setboundary(fem_mesh::FEMmesh,bdtype)

Takes in the fem_mesh and creates an array bdflag which denotes the boundary types.
1 stands for dirichlet, 2 for neumann, 3 for robin.
"""
function setboundary(node::AbstractArray,elem::AbstractArray,bdtype)
  ## Find boundary edges
  nv = size(elem,2)
  if nv == 3 # triangles
      totaledge = sort([elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]],2)
  elseif nv == 4 # quadrilateral
      totaledge = sort([elem[:,[1,2]]; elem[:,[2,3]]; elem[:,[3,4]]; elem[:,[4,1]]],2)
  end
  ne = nv # number of edges in one element
  Neall = size(totaledge,1)
  NT = size(elem,1)
  edge = unique(totaledge,1)
  totaledge = sort(totaledge,2)
  edge_matrix = sparse(round(Int,totaledge[:,1]),round(Int,totaledge[:,2]),1)
  i,j = ind2sub(size(edge_matrix),find(x->x==1,edge_matrix))
  bdedge = [i';j']'
  bdedgeidx = zeros(Int64,size(bdedge,1))
  for i = 1:size(bdedge,1)
    if VERSION < v"0.5-"
      bdedgeidx[i] = find(all(totaledge .== bdedge[i,:], 2))[1] #Find the edge in totaledge and save index
    else
      bdedgeidx[i] = find(all(totaledge .== bdedge[i,:]', 2))[1] #Find the edge in totaledge and save index
    end
  end

  bdflag = zeros(Int8,Neall)
  ## Set up boundary edges
  #nVarargin = size(varargin,2)
  #if (nVarargin==1)
  bdtype = findbdtype(bdtype)
  bdflag[bdedgeidx] = bdtype
  #end
  #=
  if (nVarargin>=2)
      for i=1:nVarargin/2
          bdtype = findbdtype(varargin{2*i-1})
          expr = varargin{2*i}
          if strcmp(expr,"all")
              bdflag(bdedgeidx) = bdtype
          else
             x = (node(allEdge(bdedgeidx,1),1) + node(allEdge(bdedgeidx,2),1))/2
             y = (node(allEdge(bdedgeidx,1),2) + node(allEdge(bdedgeidx,2),2))/2
             idx = eval(expr)
             bdflag(bdedgeidx(idx)) = bdtype
          end
      end
  end
  =#
  bdflag = reshape(bdflag,NT,ne)
  return(bdflag)
end

setboundary(fem_mesh::Mesh,bdtype) = setboundary(fem_mesh.node,fem_mesh.elem,bdtype)

function findbdtype(bdstr)
        if bdstr==:dirichlet
            bdtype = 1
        elseif bdstr==:neumann
            bdtype = 2
        elseif bdstr==:robin
            bdtype = 3
        #=
      elseif bdstr==:ABC # absorbing boundary condition for wave-type equations
            bdtype = 4
        =#
    end
    return(bdtype)
end
