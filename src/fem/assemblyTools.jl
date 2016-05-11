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
          if ~isempty(K)
            Aij = K.*Aij
          end
          A = A + sparse(elem[:,i],elem[:,j],Aij,N,N)
          if ~lumpflag
             Mij = area*((i==j)+1)/12
             M = M + sparse(elem[:,i],elem[:,j],Mij,N,N)
          end
      end
  end

  # Assemble the mass matrix by the mass lumping
  if lumpflag
      M = Diagonal(vec(accumarray([elem[:,1];elem[:,2];elem[:,3]],[area;area;area]/3,[N,1])))
  end
  return(A,M,area)
end

assemblematrix(FEMmesh;lumpflag=false,K=[]) = assemblematrix(FEMmesh.node,FEMmesh.elem,lumpflag=lumpflag,K=K)

function accumarray(subs, val, sz=(maximum(subs),))
    A = zeros(eltype(val), sz...)
    for i = 1:length(val)
        @inbounds A[subs[i]] += val[i]
    end
    A
end

function accumarray2(subs, val, fun=sum, fillval=0; sz=maximum(subs,1), issparse=false)
   counts = Dict()
   for i = 1:size(subs,1)
        counts[subs[i,:]]=[get(counts,subs[i,:],[]);val[i...]]
   end
   A = fillval*ones(sz...)
   for j = keys(counts)
        A[j...] = fun(counts[j])
   end
   issparse ? sparse(A) : A
end
