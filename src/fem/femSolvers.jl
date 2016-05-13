"""
fem_solvepoisson
"""
function fem_solvepoisson(femMesh::FEMmesh,gD::Function,f::Function,isLinear::Bool;fquadorder=3,solver="Direct",sol = [],Du = [],gN::Function = (x)->0,autodiff=true,stochastic=false,σ=(x)->0,noiseType="White")
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

  #Stochastic Part
  if stochastic
    dW = getNoise(N,node,elem,noiseType=noiseType)
    rhs(u) = quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear) + quadfbasis(σ,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear).*dW
  else #Not Stochastic
    rhs(u) = quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear)
  end
  #Solve
  if isLinear
    if solver == "Direct"
      u[freeNode]=A[freeNode,freeNode]\rhs(u)[freeNode]
    elseif solver == "CG"
      u[freeNode],ch=cg!(u[freeNode],A[freeNode,freeNode],rhs(u)[freeNode])
    elseif solver == "GMRES"
      u[freeNode],ch=gmres!(u[freeNode],A[freeNode,freeNode],rhs(u)[freeNode])
    end
    #Adjust result
    if isempty(Dirichlet) #isPureNeumann
      patchArea = accumarray(vec(elem),[area;area;area]/3, [N 1])
      uc = sum(u.*patchArea)/sum(area)
      u = u - uc   # Impose integral of u = 0
    end
  else #Nonlinear
    output = u
    function rhs!(u,resid)
      resid[freeNode]=A[freeNode,freeNode]*u[freeNode]-rhs(u)[freeNode]
    end
    nlres = nlsolve(rhs!,u,autodiff=autodiff)
    u = nlres.zero
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
"""
fem_solveheat
"""
function fem_solveheat(femMesh::FEMmesh,u0::AbstractArray,gD::Function,f::Function,isLinear::Bool;fquadorder=3,alg = "Euler",solver="LU",sol = [],Du = [],fullSave = false,saveSteps = 100,gN::FunctionOrVoid = (x,t)->0,autodiff=false,stochastic=false,σ=(x)->0,noiseType="White")
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
    tFull = Array{Float64}(round(Int64,femMesh.numIters/saveSteps))
    saveIdx = 1
    uFull[:,saveIdx] = u
  end

  #Setup for Calculations
  Minv = sparse(inv(M)) #sparse(Minv) needed until update
  if isLinear
    if alg == "Euler"
      methodType = "Explicit"
      if femMesh.μ>=0.5
        warn("Euler method chosen but μ>=.5 => Unstable. Results may be wrong.")
      end
      D = eye(N) - Δt*Minv*A
      if stochastic
        rhs(u,i,dW) = D[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode] +
                    (dW.*Minv*Δt*quadfbasis((x)->σ(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      else #Deterministic
        rhs(u,i) = D[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      end
    elseif alg == "ImplicitEuler"
      methodType = "Implicit"
      D = eye(N) + Δt*Minv*A
      lhs = D[freeNode,freeNode]
      rhs(u,i) = u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
    elseif alg == "CrankNicholson"
      methodType = "Implicit"
      Dm = eye(N) - Δt*Minv*A/2
      Dp = eye(N) + Δt*Minv*A/2
      lhs = Dp[freeNode,freeNode]
      rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
    end
  else #Nonlinear Algorithms
    if alg == "Euler"
      methodType = "Explicit"
      if femMesh.μ>=0.5
        warn("Euler method chosen but μ>=.5 => Unstable. Results may be wrong.")
      end
      D = eye(N) - Δt*Minv*A
      if stochastic
        rhs(u,i,dW) = D[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode] +
                    (dW.*Minv*Δt*quadfbasis((u,x)->σ(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      else #Deterministic
        rhs(u,i) = D[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      end
    elseif alg == "SemiImplicitEuler"
      methodType = "Implicit"
      D = eye(N) + Δt*Minv*A
      lhs = D[freeNode,freeNode]
      rhs(u,i) = u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
    elseif alg == "SemiImplicitCrankNicholson"
      methodType = "Implicit"
      Dm = eye(N) - Δt*Minv*A/2
      Dp = eye(N) + Δt*Minv*A/2
      lhs = Dp[freeNode,freeNode]
      rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
    elseif alg == "ImplicitEuler"
      methodType = "NonlinearSolve"
      function rhs!(u,output,uOld,i)
        output[freeNode] = u[freeNode] - uOld[freeNode] - Δt*Minv*A[freeNode,freeNode]*u[freeNode] - (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      end
    end
  end

  if methodType == "Implicit"
    if solver == "Cholesky"
      lhs = cholfact(lhs) # Requires positive definite, may be violated
    elseif solver == "LU"
      lhs = lufact(lhs)
    elseif solver == "QR"
      lhs = qrfact(lhs) #More stable, slower than LU
    elseif solver == "SVD"
      lhs = svdfact(lhs)
    end
  end

  if methodType == "NonlinearSolve"
    output = u
  end
  #Heat Equation Loop
  @progress for i=1:femMesh.numIters
    t = t+Δt
    if methodType == "Implicit"
      if solver == "Direct" || solver == "Cholesky" || solver == "QR" || solver == "LU" || solver == "SVD"
        u[freeNode] = lhs\rhs(u,i)
      elseif solver == "CG"
        u[freeNode],ch = cg!(u[freeNode],lhs,rhs(u,i))
      elseif solver == "GMRES"
        u[freeNode],ch = gmres!(u[freeNode],lhs,rhs(u,i))
      end
    elseif methodType == "Explicit"
      if stochastic
        dW = getNoise(N,node,elem,noiseType=noiseType)
        u[freeNode] = rhs(u,i,dW)
      else
        u[freeNode] = rhs(u,i)
      end
    elseif methodType == "NonlinearSolve"
      uOld = u
      nlres = nlsolve((u,output)->rhs!(u,output,uOld,i),u,autodiff=autodiff)
      u = nlres.zero
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
    if fullSave
      return(FEMSolution(femMesh,u,uFull,tFull))
    else
      return(FEMSolution(femMesh,u))
    end
  end
end

#fem_solveheat(femMesh::FEMmesh,u0::AbstractArray,gD::Function,f::Function,isLinear::Bool;alg = "Euler",solver="LU",sol = [],Du = [],fullSave = false,saveSteps = 100,gN::Function = (x,t)->0,autodiff=false,stochastic=false,σ=(x)->0,noiseType="White")

fem_solveheat(femMesh::FEMmesh,u0::Function,gD::Function,f::Function,isLinear::Bool;alg = "Euler",solver="LU",sol = [],Du = [],fullSave = false,saveSteps = 100,gN = (x,t)->0,autodiff=false,stochastic=false,σ=(x)->0,noiseType="White") =
fem_solveheat(femMesh::FEMmesh,u0(femMesh.node),gD::Function,f::Function,isLinear::Bool;alg = alg,solver=solver,sol = sol,Du = Du,fullSave=fullSave,saveSteps=saveSteps,gN=gN,autodiff=autodiff,stochastic=stochastic,σ=σ,noiseType=noiseType)

function fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem;alg = "Euler",solver="LU",fullSave = false,saveSteps = 100,autodiff=false,stochastic=false,σ=(x)->0,noiseType="White")
  if pdeProb.knownSol
   return(fem_solveheat(femMesh,pdeProb.u0,pdeProb.sol,pdeProb.f,pdeProb.isLinear::Bool,alg = alg,solver=solver,sol = pdeProb.sol,Du = pdeProb.Du,fullSave=fullSave,saveSteps=saveSteps,gN=pdeProb.gN,autodiff=autodiff,stochastic=pdeProb.stochastic,σ=pdeProb.σ,noiseType=pdeProb.noiseType))
  else #No Known Solution
   return(fem_solveheat(femMesh::FEMmesh,pdeProb.u0,pdeProb.gD,pdeProb.f,pdeProb.isLinear::Bool,alg = alg,solver=solver,fullSave=fullSave,saveSteps=saveSteps,gN=pdeProb.gN,autodiff=autodiff,stochastic=pdeProb.stochastic,σ=pdeProb.σ,noiseType=pdeProb.noiseType))
  end
end

function fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem;solver="Direct",fquadorder=3,autodiff=true)
  if pdeProb.knownSol
    return(fem_solvepoisson(femMesh::FEMmesh,pdeProb.gD,pdeProb.f,pdeProb.isLinear::Bool,fquadorder=fquadorder,solver=solver,sol=pdeProb.sol,Du=pdeProb.Du,gN=pdeProb.gN,autodiff=autodiff,σ=pdeProb.σ,stochastic=pdeProb.stochastic,noiseType=pdeProb.noiseType))
  else
    return(fem_solvepoisson(femMesh::FEMmesh,pdeProb.gD,pdeProb.f,pdeProb.isLinear::Bool,fquadorder=fquadorder,solver=solver,gN=pdeProb.gN,autodiff=autodiff,σ=pdeProb.σ,stochastic=pdeProb.stochastic,noiseType=pdeProb.noiseType))
  end
end

"""
quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear;gNquadorder=2)

"""
function quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear;gNquadorder=2)
  if isLinear
    bt1 = area.*(f(mid[:,:,2])+f(mid[:,:,3]))/6
    bt2 = area.*(f(mid[:,:,3])+f(mid[:,:,1]))/6
    bt3 = area.*(f(mid[:,:,1])+f(mid[:,:,2]))/6
  else
    u1 = (u[vec(elem[:,2]),:]+u[vec(elem[:,3]),:])/2
    u2 = (u[vec(elem[:,3]),:]+u[vec(elem[:,1]),:])/2
    u3 = (u[vec(elem[:,1]),:]+u[vec(elem[:,2]),:])/2
    bt1 = area.*(f(u2,mid[:,:,2])+f(u3,mid[:,:,3]))/6
    bt2 = area.*(f(u3,mid[:,:,3])+f(u1,mid[:,:,1]))/6
    bt3 = area.*(f(u1,mid[:,:,1])+f(u2,mid[:,:,2]))/6
  end
  b = vec(accumarray(vec(elem),vec([bt1;bt2;bt3])))

  if(!isempty(Dirichlet))
    uz = zeros(N)
    uz[bdNode] = gD(node[bdNode,:])
    b = b-A*uz
  end
  if(!isempty(Neumann))
    el = sqrt(float(sum((node[Neumann[:,1],:] - node[Neumann[:,2],:]).^2,2)))
    lambdagN,weightgN = quadpts1(gNquadorder)
    phigN = lambdagN                # linear bases
    nQuadgN = size(lambdagN,1)
    ge = zeros(size(Neumann,1),2)
    for pp = 1:nQuadgN
        # quadrature points in the x-y coordinate
        ppxy = lambdagN[pp,1]*node[Neumann[:,1],:] +
               lambdagN[pp,2]*node[Neumann[:,2],:]
        gNp = gN(ppxy)
        for igN = 1:2
            ge[:,igN] = ge[:,igN] + weightgN[pp]*phigN[pp,igN]*gNp
        end
    end
    ge = ge.*repmat(el,1,2)
    b = b + accumarray(vec(Neumann), vec(ge),[N,1])
  end
  return(b)
end
