macro def(name, definition)
    return quote
        macro $name()
            esc($(Expr(:quote, definition)))
        end
    end
end

## Unused other versions of Functions

#=
"""
quadfbasis2(f,gD,A,node,elem,lambda,phi,weight,N,NT,area,bdnode)
Slightly slower than quadfbasis, easier to extend to higher order quadrature
"""
function quadfbasis2(f,gD,A,node,elem,lambda,phi,weight,N,NT,area,bdnode)
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
  if(!isempty(dirichlet))
    uz = zeros(N)
    uz[bdnode] = gD(node[bdnode,:])
    b = b-A*uz
    if(!isempty(neumann))
      Nve = node[neumann[:,1],:] - node[neumann[:,2],:]
      edgeLength = sqrt(sum(Nve.^2,2))
      mid = (node[neumann[:,1],:] + node[neumann[:,2],:])/2
      b = b + accumarray(int([vec(neumann),ones(2*size(neumann,1),1)]), repmat(edgeLength.*g_N(mid)/2,2,1),[N,1])
    end
  else #Pure neumann
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

"""
Splats keys from a dict into variables

```
@materialize a, b, c = dict
```

"""
macro materialize(dict_splat)
    keynames, dict = dict_splat.args
    keynames = isa(keynames, Symbol) ? [keynames] : keynames.args
    dict_instance = gensym()
    kd = [:($key = $dict_instance[$(Expr(:quote, key))]) for key in keynames]
    kdblock = Expr(:block, kd...)
    expr = quote
        $dict_instance = $dict # handle if dict is not a variable but an expression
        $kdblock
    end
    esc(expr)
end

function shapeResult(res)
  #Input is a vector of tuples
  #Output the "columns" as vectors
  out = cell(length(res[1]))
  for i = 1:length(res[1])
    out[i] = Vector{typeof(res[1][i])}(length(res))
  end
  for i = 1:length(res), j = 1:length(out)
    out[j][i] = res[i][j]
  end
  return tuple(out...)
end

function monteCarlo(func::Function,args...)
  res = pmap(func,args...)
  if length(res[1])==1
    return res
  else
    return shapeResult(res)
  end
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
accumarray(subs, val, sz=(maximum(subs),))

See MATLAB's documentation for more details.
"""
function accumarray(subs, val, sz=(maximum(subs),))
    A = zeros(eltype(val), sz...)
    for i = 1:length(val)
        @inbounds A[subs[i]] += val[i]
    end
    A
end

"""
Slower than accumarray but more functionality
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
