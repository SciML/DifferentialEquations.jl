"""
fem_solvepoisson
"""
function fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem;solver::AbstractString="Direct",fquadorder::Int64=3,autodiff::Bool=true)
  #Assemble Matrices
  A,M,area = assemblematrix(femMesh,lumpflag=true)

  #Unroll some important constants
  @unpack femMesh: Δt,bdNode,node,elem,N,NT,freeNode,Dirichlet,Neumann
  @unpack pdeProb: f,Du,f,gD,gN,sol,knownSol,isLinear,σ,stochastic,noiseType

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
  if knownSol # True solution exists
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

`fem_solveheat(femMesh::FEMmesh,u0::AbstractArray,gD::Function,f::Function,isLinear::Bool)`

`fem_solveheat(femMesh::FEMmesh,u0::Function,gD::Function,f::Function,isLinear::Bool)`

`fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem)`

Takes in a definition for the heat equation ``u_t = Δu + f`` on a finite element
mesh with initial condtion u0 and returns the solution.
"""
function fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem;alg::AbstractString = "Euler",fquadorder::Int64=3,
  solver::AbstractString="LU",fullSave::Bool = false,saveSteps::Int64 = 100,autodiff::Bool=false)
  #Assemble Matrices
  A,M,area = assemblematrix(femMesh,lumpflag=true)

  #Unroll some important constants
  @unpack femMesh: Δt,bdNode,node,elem,N,NT,freeNode,Dirichlet,Neumann
  @unpack pdeProb: f,u0,Du,f,gD,gN,sol,knownSol,isLinear,σ,stochastic,noiseType

  #Note if Atom is loaded for progress
  atomLoaded = checkIfLoaded("Atom")

  #Set Initial
  u = u0(node)
  t = 0
  #Setup f quadraturef
  mid = Array{Float64}(size(node[vec(elem[:,2]),:])...,3)
  mid[:,:,1] = (node[vec(elem[:,2]),:]+node[vec(elem[:,3]),:])/2
  mid[:,:,2] = (node[vec(elem[:,3]),:]+node[vec(elem[:,1]),:])/2
  mid[:,:,3] = (node[vec(elem[:,1]),:]+node[vec(elem[:,2]),:])/2

  #Setup animation

  if fullSave
    uFull = Array{Float64}(length(u),ceil(Int64,femMesh.numIters/saveSteps))
    tFull = Array{Float64}(ceil(Int64,femMesh.numIters/saveSteps))
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
      if stochastic
        rhs(u,i,dW) = u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode] +
                    (dW.*Minv*Δt*quadfbasis((x)->σ(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      else #Deterministic
        rhs(u,i) = u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      end
    elseif alg == "CrankNicholson"
      methodType = "Implicit"
      Dm = eye(N) - Δt*Minv*A/2
      Dp = eye(N) + Δt*Minv*A/2
      lhs = Dp[freeNode,freeNode]
      if stochastic
        rhs(u,i,dW) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode] +
                    (dW.*Minv*Δt*quadfbasis((x)->σ(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      else #Deterministic
        rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((x)->f(x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      end
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
      if stochastic
        rhs(u,i,dW) = u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode] +
                    (dW.*Minv*Δt*quadfbasis((u,x)->σ(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      else #Deterministic
        rhs(u,i) = u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      end
    elseif alg == "SemiImplicitCrankNicholson"
      methodType = "Implicit"
      Dm = eye(N) - Δt*Minv*A/2
      Dp = eye(N) + Δt*Minv*A/2
      lhs = Dp[freeNode,freeNode]
      if stochastic
        rhs(u,i,dW) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode] +
                    (dW.*Minv*Δt*quadfbasis((u,x)->σ(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      else #Deterministic
        rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
      end
    elseif alg == "ImplicitEuler"
      methodType = "NonlinearSolve"
      if stochastic
        function rhs!(u,output,dW,uOld,i)
          output[freeNode] = u[freeNode] - uOld[freeNode] - Δt*Minv*A[freeNode,freeNode]*u[freeNode] - (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode] -(dW.*Minv*Δt*quadfbasis((u,x)->σ(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                      A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
        end
      else #Deterministic
        function rhs!(u,output,uOld,i)
          output[freeNode] = u[freeNode] - uOld[freeNode] - Δt*Minv*A[freeNode,freeNode]*u[freeNode] - (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear))[freeNode]
        end
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
  for i=1:femMesh.numIters
    t += Δt
    if methodType == "Implicit"
      if stochastic
        dW = getNoise(N,node,elem,noiseType=noiseType)
        if solver == "Direct" || solver == "Cholesky" || solver == "QR" || solver == "LU" || solver == "SVD"
          u[freeNode] = lhs\rhs(u,i,dW)
        elseif solver == "CG"
          u[freeNode],ch = cg!(u[freeNode],lhs,rhs(u,i,dW))
        elseif solver == "GMRES"
          u[freeNode],ch = gmres!(u[freeNode],lhs,rhs(u,i,dW))
        end
      else #Deterministic
        if solver == "Direct" || solver == "Cholesky" || solver == "QR" || solver == "LU" || solver == "SVD"
          u[freeNode] = lhs\rhs(u,i)
        elseif solver == "CG"
          u[freeNode],ch = cg!(u[freeNode],lhs,rhs(u,i))
        elseif solver == "GMRES"
          u[freeNode],ch = gmres!(u[freeNode],lhs,rhs(u,i))
        end
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
      if stochastic
        dW = getNoise(N,node,elem,noiseType=noiseType)
        nlres = nlsolve((u,output)->rhs!(u,output,dW,uOld,i),u,autodiff=autodiff)
      else
        nlres = nlsolve((u,output)->rhs!(u,output,uOld,i),u,autodiff=autodiff)
      end
      u = nlres.zero
    end
    u[bdNode] = gD(node,i*Δt)[bdNode]
    if fullSave && i%saveSteps==0
      saveIdx+=1
      uFull[:,saveIdx] = u
      tFull[saveIdx] = t
    end
    atomLoaded ? progress(i/femMesh.numIters) : nothing
  end
  if knownSol #True Solution exists
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
