"""
assemblematrix(node,elem;lumpflag=false,K=[])

Assembles the stiffness matrix A as an approximation to Δ
on the finite element mesh (node,elem). Also generates the
mass matrix M. If lumpflag=true, then the mass matrix is lumped
resulting in a diagonal mass matrix. Specify a diffusion constant
along the nodes via K.

### Returns
A = Stiffness Matrix
M = Mass Matrix
area = A vector of the calculated areas for each element.
"""
function assemblematrix(node,elem;lumpflag=false,K=[])
  ## ASSEMBLEMATRIX matrix for diffusion and reaction

  # Parameters
  N = size(node,1)
  A = spzeros(N,N)
  M = spzeros(N,N)

  # 3-D case
  if (size(node,2) == 3) && (size(elem,2) == 4) # 3-D
      A,M,area = assemblematrix3(node,elem,lumpflag) #Not Implemented
      return
  end

  ve = Array{Float64}(size(node[elem[:,3],:])...,3)
  ## Compute vedge, edge as a vector, and area of each element
  ve[:,:,1] = node[elem[:,3],:]-node[elem[:,2],:]
  ve[:,:,2] = node[elem[:,1],:]-node[elem[:,3],:]
  ve[:,:,3] = node[elem[:,2],:]-node[elem[:,1],:]
  area = 0.5*abs(-ve[:,1,3].*ve[:,2,2]+ve[:,2,3].*ve[:,1,2])

  # Assemble stiffness matrix
  for i = 1:3
      for j = 1:3
          Aij = (ve[:,1,i].*ve[:,1,j]+ve[:,2,i].*ve[:,2,j])./(4*area)
          if !isempty(K)
            Aij = K.*Aij
          end
          A = A + sparse(elem[:,i],elem[:,j],Aij,N,N)
          if !lumpflag
             Mij = area*((i==j)+1)/12
             M = M + sparse(elem[:,i],elem[:,j],Mij,N,N)
          end
      end
  end

  # Assemble the mass matrix by the mass lumping
  if lumpflag
      M = Diagonal(vec(Matlab.accumarray([elem[:,1];elem[:,2];elem[:,3]],[area;area;area]/3,[N,1])))
  end
  return(A,M,area)
end

"""
assemblematrix(FEMmesh::FEMmesh;lumpflag=false,K=[])

Assembles the stiffness matrix A as an approximation to Δ
on the finite element mesh (node,elem). Also generates the
mass matrix M. If lumpflag=true, then the mass matrix is lumped
resulting in a diagonal mass matrix. Specify a diffusion constant
along the nodes via K.

### Returns
A = Stiffness Matrix
M = Mass Matrix
area = A vector of the calculated areas for each element.
"""
assemblematrix(FEMmesh::FEMmesh;lumpflag=false,K=[]) = assemblematrix(FEMmesh.node,FEMmesh.elem,lumpflag=lumpflag,K=K)
