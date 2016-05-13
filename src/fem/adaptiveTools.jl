"""
estimaterecovery(node,elem,u)
"""
function estimaterecovery(node,elem,u)
  #Computes the Δ error estimate η
  Dlambda,area = gradbasis(node,elem)
  Du = gradu(node,elem,u,Dlambda)
  Du = recovery(node,elem,Du,area)
  DDu[:,1:2] = gradu(node,elem,Du[:,1],Dlambda)
  DDu[:,3:4] = gradu(node,elem,Du[:,2],Dlambda)
  η = area.*sum(abs(DDu),2)
  return(η,Du)
end
"""
recovery(node,elem,Du,area)
"""
function recovery(node,elem,Du,area)
  #Promotes a piecewise constant function to piecewise linear
  N = size(node,1)
  dudxArea = area.*Du[:,1]
  dudyArea = area.*Du[:,2]
  patchArea = accumarray(vec(elem),[area;area;area], [N 1])
  dudxArea = accumarray(vec(elem),[dudxArea;dudxArea;dudxArea],[N 1])
  dudyArea = accumarray(vec(elem),[dudyArea;dudyArea;dudyArea],[N 1])
  dudx = dudxArea./patchArea
  dudy = dudyArea./patchArea
  RDu = [dudx dudy]
  return(RDu)
end

"""
mark(elem,eta,theta;method="L2")
"""
function mark(elem,eta,theta;method="L2")
  NT = size(elem,1); isMark = false(NT,1)
  if method == "Max"
    isMark[eta>theta*max(eta)]=1
  elseif method == "COARSEN"
    isMark[eta<theta*max(eta)]=1
  elseif method == "L2"
    sortedEta,idx = sort(eta.^2,rev=true)
    x = cumsum(sortedEta)
    isMark[idx[x < theta* x[NT]]] = 1
    isMark[idx[1]] = 1
  end
  markedElem = convert(Int64,find(isMark==true))
  return(markedElem)
end

"""
bisect(node,elem;markedElem= (1:size(elem,1))',bdFlag)
"""
function bisect(node,elem;markedElem=1:size(elem,1),bdFlag=[])
  # BISECT bisect a 2-D triangulation.

  # Set up
  HB = []; tree = []
  if isempty(markedElem) return(node,elem,bdFlag,HB,tree) end
  if markedElem == "all" markedElem = 1:size(elem,1) end
  if isa(markedElem,Array{Bool}) markedElem = convert(Int64,markedElem) end

  # Construct auxiliary data structure
  T = auxstructure(elem)
  neighbor = T.neighbor; elem2edge = T.elem2edge; edge = T.edge
  #clear T;
  #[neighbor,elem2edge,edge] = auxstructurec(int32(elem));
  N = size(node,1); NT = size(elem,1); NE = size(edge,1)

  # Add new nodes
  isCutEdge = false(NE,1)
  while sum(markedElem)>0
      #isCutEdge[elem2edge(markedElem,1)] = true
      refineNeighbor = neighbor(markedElem,1)
      markedElem = refineNeighbor(~isCutEdge(elem2edge(refineNeighbor,1)))
  end
  edge2newNode = zeros(Int64,NE,1)
  edge2newNode[isCutEdge] = N+1:N+sum(isCutEdge)
  HB = zeros(Int64,sum(isCutEdge),3)
  HB[:,1] = edge2newNode(isCutEdge)
  HB[:,[2 3]] = edge(isCutEdge,[1 2])
  node[HB[:,1],:] = [node[HB[:,2],:] + node[HB[:,3],:]]/2;

  # Refine marked elements
  Nb = 0; tree = zeros(Int64,3*NT,3)
  for k = 1:2
      t = find(edge2newNode(elem2edge(:,1))>0)
      newNT = length(t)
      if (newNT == 0) break; end
      L = t; R = NT+1:NT+newNT
      p1 = elem(t,1); p2 = elem(t,2); p3 = elem(t,3)
      p4 = edge2newNode(elem2edge(t,1))
      elem[L,:] = [p4 p1 p2]
      elem[R,:] = [p4 p3 p1]
    if nargin==4 && ~isempty(bdFlag) # Refine boundary edges
     		#bdFlag[R [1 3]] = bdFlag[t [2 1]]
     		#bdFlag[L [1 2]] = bdFlag[t [3 1]]
          bdFlag[L,3] = 0
      else
          bdFlag = []
  	end
      tree[Nb+1:Nb+newNT,1] = L
      tree[Nb+1:Nb+newNT,2] = L
      tree[Nb+1:Nb+newNT,3] = R
      elem2edge[L,1] = elem2edge[t,3]
      elem2edge[R,1] = elem2edge[t,2]
      NT = NT + newNT; Nb = Nb + newNT
  end
  tree = tree[1:Nb,:]
  return(node,elem,bdFlag,HB,tree)
end
