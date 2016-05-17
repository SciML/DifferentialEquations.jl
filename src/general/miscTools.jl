"""
`modulechildren(m::Module)`

Returns the modules in m
"""
modulechildren(m::Module) = filter(x->isa(x, Module), map(x->m.(x), names(m, true)))

"""
`checkIfLoaded(pkg::AbstractString)`

Returns true if module "pkg" is defined in Main, otherwise false.
"""
checkIfLoaded(pkg::AbstractString)= maximum(map(string,modulechildren(Main)).==pkg)

## Unused other versions of Functions

#=
"""
quadfbasis2(f,gD,A,node,elem,lambda,phi,weight,N,NT,area,bdNode)
Slightly slower than quadfbasis, easier to extend to higher order quadrature
"""
function quadfbasis2(f,gD,A,node,elem,lambda,phi,weight,N,NT,area,bdNode)
  nQuad = size(lambda,1)
  bt = zeros(NT,3)
  for p = 1:nQuad
      pxy = lambda[p,1]*node[elem[:,1],:] +
        lambda[p,2]*node[elem[:,2],:] +
        lambda[p,3]*node[elem[:,3],:]
      fp = f(pxy)
      for i = 1:3
          bt[:,i] = bt[:,i] + weight[p]*phi[p,i]*fp
      end
  end
  bt = bt.*repmat(area,1,3)
  b = vec(accumarray(vec(elem),vec(bt),[N 1]))
  if(!isempty(Dirichlet))
    uz = zeros(N)
    uz[bdNode] = gD(node[bdNode,:])
    b = b-A*uz
    if(!isempty(Neumann))
      Nve = node[Neumann[:,1],:] - node[Neumann[:,2],:]
      edgeLength = sqrt(sum(Nve.^2,2))
      mid = (node[Neumann[:,1],:] + node[Neumann[:,2],:])/2
      b = b + accumarray(int([vec(Neumann),ones(2*size(Neumann,1),1)]), repmat(edgeLength.*g_N(mid)/2,2,1),[N,1])
    end
  else #Pure Neumann
    b = b-mean(b) #Compatibility condition: sum(b)=0
    b[1] = 0 #Fix 1 point
  end
  return(b)
end

"""
CG2(u,A,b;tol=1e-6)

Needs to be tested. Could be faster than CG.
"""
function CG2(u,A,b;tol=1e-6)
  tol = tol*norm(b)
  k = 1
  r = b - A*u
  p = r
  r2 = dot(r,r)
  while sqrt(r2) >= tol && k<length(b)
    Ap = A*p
    alpha = r2/dot(p,Ap)
    u = u + alpha*p
    r = r - alpha*Ap
    r2old = r2
    r2 = dot(r,r)
    beta = r2/r2old
    p = r + beta*p
    k = k + 1
  end
  return u,k
end
=#

#=
"""
Slower than accumarray
"""
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
=#
