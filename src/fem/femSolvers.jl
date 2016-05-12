
function fem_solvepoisson(femMesh::FEMmesh,gD::Function,gN::Function,f::Function;fquadorder=3,solver="Direct",sol = [],Du = [])
  #Assemble Matrices
  A,M,area = assemblematrix(femMesh,lumpflag=true)

  #Unroll some important constants
  bdNode = femMesh.bdNode
  node = femMesh.node
  elem = femMesh.elem
  N = femMesh.N
  NT = femMesh.NT
  freeNode = femMesh.freeNode
  Dirichlet = femMesh.Dirichlet
  Neumann = femMesh.Neumann

  #Setup f quadrature
  mid = Array{Float64}(size(node[vec(elem[:,2]),:])...,3)
  mid[:,:,1] = (node[vec(elem[:,2]),:]+node[vec(elem[:,3]),:])/2
  mid[:,:,2] = (node[vec(elem[:,3]),:]+node[vec(elem[:,1]),:])/2
  mid[:,:,3] = (node[vec(elem[:,1]),:]+node[vec(elem[:,2]),:])/2

  #Setup u
  u = zeros(N)
  u[bdNode] = gD(node[bdNode,:])

  #Solve
  if solver == "Direct"
    u[freeNode]=A[freeNode,freeNode]\quadfbasis(f,gD,gN,A,node,elem,area,bdNode,mid,N,Dirichlet,Neumann)[freeNode]
  elseif solver == "CG"
    u[freeNode],ch=cg!(u[freeNode],A[freeNode,freeNode],quadfbasis(f,gD,gN,A,node,elem,area,bdNode,mid,N,Dirichlet,Neumann)[freeNode])
  elseif solver == "GMRES"
    u[freeNode],ch=gmres!(u[freeNode],A[freeNode,freeNode],quadfbasis(f,gD,gN,A,node,elem,area,bdNode,mid,N,Dirichlet,Neumann)[freeNode])
  end
  #Adjust result
  if isempty(Dirichlet) #isPureNeumann
    patchArea = accumarray(vec(elem),[area;area;area]/3, [N 1])
    uc = sum(u.*patchArea)/sum(area)
    u = u - uc   # Impose integral of u = 0
  end
  #Return
  if typeof(sol)==Function # True solution exists
    return(FEMSolution(femMesh,u,sol(node),sol,Du))
  else #No true solution
    return(FEMSolution(femMesh,u))
  end
end

## Evolution Equation Solvers
#Note
#rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode] + Δt*f(node,(i-.5)*Δt)[freeNode] #Nodel interpolation 1st order

function fem_solveheat(femMesh::FEMmesh,u0::AbstractArray,gD::Function,f::Function;fquadorder=3,alg = "Euler",solver="CG",sol = [],Du = [],fullSave = false,saveSteps = 100)
  #Assemble Matrices
  A,M,area = assemblematrix(femMesh,lumpflag=true)

  #Unroll some important constants
  Δt  = femMesh.Δt
  bdNode = femMesh.bdNode
  node = femMesh.node
  elem = femMesh.elem
  N = femMesh.N
  NT = femMesh.NT
  freeNode = femMesh.freeNode
  Dirichlet = femMesh.Dirichlet
  Neumann = femMesh.Neumann

  #Set Initial
  u = u0
  t = 0
  #Setup f quadrature
  #=
  #For quadfbasis2
  lambda,weight = quadpts(fquadorder)
  phi = lambda
  =#
  mid = Array{Float64}(size(node[vec(elem[:,2]),:])...,3)
  mid[:,:,1] = (node[vec(elem[:,2]),:]+node[vec(elem[:,3]),:])/2
  mid[:,:,2] = (node[vec(elem[:,3]),:]+node[vec(elem[:,1]),:])/2
  mid[:,:,3] = (node[vec(elem[:,1]),:]+node[vec(elem[:,2]),:])/2

  #Setup animation

  if fullSave
    uFull = Array{Float64}(length(u),round(Int64,femMesh.numIters/saveSteps))
    saveIdx = 1
    uFull[:,aniI] = u
  end

  #Setup for Calculations
  Minv = sparse(inv(M)) #sparse(Minv) needed until update
  if alg == "Euler"
    implicitMethod = false
    if femMesh.μ>=0.5
      warn("Euler method chosen but μ>=.5 => Unstable. Results may be wrong.")
    end
    D = eye(N) - Δt*Minv*A
    rhs(u,i) = D[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),A,node,elem,area,bdNode,mid,N,Dirichlet,Neumann))[freeNode]
  elseif alg == "ImplicitEuler"
    implicitMethod = true
    D = eye(N) + Δt*Minv*A
    lhs = D[freeNode,freeNode]
    rhs(u,i) = u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,node,elem,area,bdNode,mid,N,Dirichlet,Neumann))[freeNode]
  elseif alg == "CrankNicholson"
    implicitMethod = true
    Dm = eye(N) - Δt*Minv*A/2
    Dp = eye(N) + Δt*Minv*A/2
    lhs = Dp[freeNode,freeNode]
    rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),A,node,elem,area,bdNode,mid,N,Dirichlet,Neumann))[freeNode]
  end

  if solver == "Cholesky" && implicitMethod
    lhs = cholfact(lhs) # Requires positive definite, may be violated
  elseif solver == "LU"
    lhs = lufact(lhs)
  elseif solver == "QR"
    lhs = qrfact(lhs) #More stable, slower than LU
  elseif solver == "SVD"
    lhs = svdfact(lhs)
  end

  #Heat Equation Loop
  @progress for i=1:femMesh.numIters
    t = t+Δt
    if implicitMethod
      if solver == "Direct" || solver == "Cholesky" || solver == "QR" || solver == "LU" || solver == "SVD"
        u[freeNode] = lhs\rhs(u,i)
      elseif solver == "CG"
        u[freeNode],ch = cg!(u[freeNode],lhs,rhs(u,i))
      elseif solver == "GMRES"
        u[freeNode],ch = gmres!(u[freeNode],lhs,rhs(u,i))
      end
    else #explicitMethod
      u[freeNode] = rhs(u,i)
    end
    u[bdNode] = gD(node,i*Δt)[bdNode]
    if fullSave && i%saveSteps==0
      saveIdx+=1
      uFull[:,saveIdx] = u
      tFull[saveIdx] = t
    end
  end
  if typeof(sol)==Function #True Solution exists
    if fullSave
      return(FEMSolution(femMesh,u,sol(node,femMesh.T),sol,Du,uFull,tFull))
    else
      return(FEMSolution(femMesh,u,sol(node,femMesh.T),sol,Du))
    end
  else #No true solution
    return(FEMSolution(femMesh,u))
  end
end

fem_solveheat(femMesh::FEMmesh,u0::Function,gD::Function,f::Function;alg = "Euler",solver="CG",sol = [],Du = [],fullSave = false,saveSteps = 100) = fem_solveheat(femMesh::FEMmesh,u0(femMesh.node),gD::Function,f::Function;alg = alg,solver=solver,sol = sol,Du = Du,animate=animate,animateSteps=animateSteps)

function fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem;alg = "Euler",solver="CG",fullSave = false,saveSteps = 100)
  if pdeProb.knownSol
   return(fem_solveheat(femMesh,pdeProb.u0,pdeProb.sol,pdeProb.f,alg = alg,solver=solver,sol = pdeProb.sol,Du = pdeProb.Du,fullSave=fullSave,saveSteps=saveSteps))
  else #No Known Solution
   return(fem_solveheat(femMesh::FEMmesh,pdeProb.u0,pdeProb.gD,pdeProb.f,alg = alg,solver=solver,fullSave=fullSave,saveSteps=saveSteps))
  end
end

function fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem;solver="Direct",fquadorder=3)
  if pdeProb.knownSol
    return(fem_solvepoisson(femMesh::FEMmesh,pdeProb.gD,pdeProb.gN,pdeProb.f,fquadorder=fquadorder,solver=solver,sol=pdeProb.sol,Du=pdeProb.Du))
  else
    return(fem_solvepoisson(femMesh::FEMmesh,pdeProb.gD,pdeProb.gN,pdeProb.f,fquadorder=fquadorder,solver=solver))
  end
end

function quadfbasis(f,gD,gN,A,node,elem,area,bdNode,mid,N,Dirichlet,Neumann)
  bt1 = area.*(f(mid[:,:,2])+f(mid[:,:,3]))/6
  bt2 = area.*(f(mid[:,:,3])+f(mid[:,:,1]))/6
  bt3 = area.*(f(mid[:,:,1])+f(mid[:,:,2]))/6
  b = vec(accumarray(vec(elem),vec([bt1;bt2;bt3])))

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
  else #Pure Neumann 0
    b = b-mean(b) #Compatibility condition: sum(b)=0
    b[1] = 0 #Fix 1 point
  end
  return(b)
end

function quadfbasis2(f,gD,A,node,elem,lambda,phi,weight,N,NT,area,bdNode)
  #for Crank-Nicholson
  #rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis2((x)->f(x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),A,node,elem,lambda,phi,weight,N,NT,area,bdNode))[freeNode]
  #Slightly slower, easier to extend to higher order quadrature
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
