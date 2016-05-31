"""
## Finite Element Poisson Equation Solver

solve(femMesh::FEMmesh,pdeProb::PoissonProblem)

Takes in a definition for the heat equation ``-Δu = f`` on `femMesh` with
functions as defined in `pdeProb`. If `σ` is specified in `pdeProb`, then this
solves the stochastic Poisson equation ``-Δu = f + σdW``.

### Keyword Arguments

* `solver` = Linear solver algorithm. This is the algorithm which is chosen for solving
the implicit equation `Ax=b`. The default is `LU`. The choices are:
  * `Direct` = Solves `Ax=b` using `\`
  * `CG` = Conjugate-Gradient. Best when the space is very large and ``I ± ΔtM⁻¹A`` is positive definite.
  * `GMRES` = GMRES. Best when the space is very large and ``I ± ΔtM⁻¹A`` is not positive definite.
* `saveSteps` = If `fullSave=true`, then this is the number of steps between the saves.
* `autodiff` = Whether or not autodifferentiation (as provided by AutoDiff.jl) is used
for the nonlinear solving. By default autodiff is false.
"""
function solve(femMesh::FEMmesh,pdeProb::PoissonProblem;solver::String="Direct",autodiff::Bool=false,method=:trust_region,show_trace=false,iterations=1000)
  #Assemble Matrices
  A,M,area = assemblematrix(femMesh,lumpflag=true)

  #Unroll some important constants
  @unpack femMesh: Δt,bdNode,node,elem,N,NT,freeNode,Dirichlet,Neumann
  @unpack pdeProb: f,Du,f,gD,gN,sol,knownSol,isLinear,u₀,numVars,σ,stochastic,noiseType, D

  #Setup f quadrature
  mid = Array{Float64}(size(node[vec(elem[:,2]),:])...,3)
  mid[:,:,1] = (node[vec(elem[:,2]),:]+node[vec(elem[:,3]),:])/2
  mid[:,:,2] = (node[vec(elem[:,3]),:]+node[vec(elem[:,1]),:])/2
  mid[:,:,3] = (node[vec(elem[:,1]),:]+node[vec(elem[:,2]),:])/2

  #Setup u
  u = u₀(node)
  if numVars==0
    numVars = size(u,2)
    #pdeProb.numVars = numVars #Mutate problem to be correct.
    if gD == nothing
      gD=(x)->zeros(size(x,1),numVars)
    end
    if gN == nothing
      gN=(x)->zeros(size(x,1),numVars)
    end
    if D == nothing
      if numVars == 1
        D = 1.0
      else
        D = ones(1,numVars)
      end
    end
  end
  u[bdNode] = gD(node[bdNode,:])

  #Stochastic Part
  if stochastic
    rands = getNoise(u,node,elem,noiseType=noiseType)
    dW = next(rands)
    rhs(u) = quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars) + quadfbasis(σ,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars).*dW
  else #Not Stochastic
    rhs(u) = quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars)
  end
  #Solve
  if isLinear
    if solver == "Direct"
      u[freeNode,:]=D.*(A[freeNode,freeNode]\rhs(u)[freeNode])
    elseif solver == "CG"
      u[freeNode,:],ch=cg!(u[freeNode,:],A[freeNode,freeNode],rhs(u)[freeNode]) # Needs diffusion constant
    elseif solver == "GMRES"
      u[freeNode,:],ch=gmres!(u[freeNode,:],A[freeNode,freeNode],rhs(u)[freeNode]) # Needs diffusion constants
    end
    #Adjust result
    if isempty(Dirichlet) #isPureNeumann
      patchArea = accumarray(vec(elem),[area;area;area]/3, [N 1])
      uc = sum(u.*patchArea)/sum(area)
      u = u - uc   # Impose integral of u = 0
    end
  else #Nonlinear
    function rhs!(u,resid)
      u = reshape(u,N,numVars)
      resid = reshape(resid,N,numVars)
      resid[freeNode,:]=D.*(A[freeNode,freeNode]*u[freeNode,:])-rhs(u)[freeNode,:]
      u = vec(u)
      resid = vec(resid)
    end
    u = vec(u)
    nlres = nlsolve(rhs!,u,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
    u = nlres.zero
    if numVars > 1
      u = reshape(u,N,numVars)
    end
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
#rhs(u,i) = Dm[freeNode,freeNode]*u[freeNode,:] + Δt*f(node,(i-.5)*Δt)[freeNode] #Nodel interpolation 1st 𝒪
"""
## Finite Element Heat Equation Solver

`solve(femMesh::FEMmesh,pdeProb::HeatProblem)`

Takes in a definition for the heat equation ``u_t = Δu + f`` on `femMesh` with
functions as defined in `pdeProb`. If `σ` is specified in `pdeProb`, then this
solves the stochastic heat equation ``u_t = Δu + f + σdW_t``.

### Keyword Arguments

* `alg` = Solution algorithm. Default is Euler. The choices are:
  * Linear
    * Euler (Explicit)
    * Implicit Euler (Implicit)
    * Crank-Nicholson (Implicit)
  * Nonlinear
    * Euler (Explicit)
    * Implicit Euler (Nonlinear Solve)
    * Crank-Nicholson (Nonlinear Solve)
    * Semi-Implicit Euler (Implicit)
    * Semi-Implicit Crank-Nicholson (Implicit)

Explicit algorithms only require solving matrix multiplications `Au`. Implicit algorithms
require solving the linear equation `Ax=b` where `x` is the unknown. Nonlinear Solve algorithms
require solving the nonlinear equation f(x)=0 using methods like Newton's method and is
provided by NLSolve.jl. Explicit algorithms have the least stability and should
be used either small Δt and non-stiff equations. The implicit algorithms have better stability,
but for nonlinear equations require costly nonlinear solves in order to be solved exactly.
The semi-implicit algorithms discretize with part of the equation implicit and another
part explicit in order to allow for the algorithm to not require a nonlinear solve, but
at the cost of some stability (though still vastly better at stability than explicit algorithms).

* `solver` = Linear solver algorithm. This is the algorithm which is chosen for solving
the implicit equation `Ax=b`. The default is `LU`. The choices are:
  * `Direct` = Solves using `\` (no factorization). Not recommended.
  * `Cholesky` = Cholsky decomposition. Only stable of ``I ± ΔtM⁻¹A`` is positive definite.
    This means that this works best when Δt is small. When applicable, this is the fastest.
  * `LU` = LU-Decomposition. A good mix between fast and stable.
  * `QR` = QR-Decomposition. Less numerical roundoff error than `LU`, but slightly slower.
  * `SVD` = SVD-Decomposition. By far the slowest, but the most robust to roundoff error.
  * `CG` = Conjugate-Gradient. Best when the space is very large and ``I ± ΔtM⁻¹A`` is positive definite.
  * `GMRES` = GMRES. Best when the space is very large and ``I ± ΔtM⁻¹A`` is not positive definite.
* `fullSave` = Makes the algorithm save the output at every `saveSteps` timesteps.
By default fullSave is false.
* `saveSteps` = If `fullSave=true`, then this is the number of steps between the saves.
* `autodiff` = Whether or not autodifferentiation (as provided by AutoDiff.jl) is used
for the nonlinear solving. By default autodiff is false.
"""
function solve(femMesh::FEMmesh,pdeProb::HeatProblem;alg::String = "Euler",
  solver::AbstractString="LU",fullSave::Bool = false,saveSteps::Int = 100,
  autodiff::Bool=false,method=:trust_region,show_trace=false,iterations=1000)
  #Assemble Matrices
  A,M,area = assemblematrix(femMesh,lumpflag=true)

  #Unroll some important constants
  @unpack femMesh: Δt,bdNode,node,elem,N,NT,freeNode,Dirichlet,Neumann
  @unpack pdeProb: f,u₀,Du,gD,gN,sol,knownSol,isLinear,numVars,σ,stochastic,noiseType,D

  #Note if Atom is loaded for progress
  atomLoaded = isimported("Atom")

  #Set Initial
  u = u₀(node)
  if numVars==0
    numVars = size(u,2)
    #pdeProb.numVars = numVars #Mutate problem to be correct.
    if gD == nothing
      gD=(x,t)->zeros(size(x,1),numVars)
    end
    if gN == nothing
      gN=(x,t)->zeros(size(x,1),numVars)
    end
    if D == nothing
      if numVars == 1
        D = 1.0
      else
        D = ones(1,numVars)
      end
    end
  end
  t = 0
  √Δt = sqrt(Δt)
  #Setup f quadraturef
  mid = Array{Float64}(size(node[vec(elem[:,2]),:])...,3)
  mid[:,:,1] = (node[vec(elem[:,2]),:]+node[vec(elem[:,3]),:])/2
  mid[:,:,2] = (node[vec(elem[:,3]),:]+node[vec(elem[:,1]),:])/2
  mid[:,:,3] = (node[vec(elem[:,1]),:]+node[vec(elem[:,2]),:])/2

  #Setup animation

  if fullSave
    uFull = GrowableArray(u)
    tFull = Float64[t]
  end

  #Setup for Calculations
  Minv = sparse(inv(M)) #sparse(Minv) needed until update
  if isLinear
    if alg == "Euler"
      methodType = "Explicit"
      if femMesh.μ>=0.5
        warn("Euler method chosen but μ>=.5 => Unstable. Results may be wrong.")
      end
      K = eye(N) - Δt*Minv*D*A #D okay since numVar = 1 for linear
      if stochastic
        rhs(u,i,dW) = K[freeNode,freeNode]*u[freeNode,:] + (Minv*Δt*quadfbasis((x)->f(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:] +
                    (√Δt.*dW.*Minv*quadfbasis((x)->σ(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      else #Deterministic
        rhs(u,i) = K[freeNode,freeNode]*u[freeNode,:] + (Minv*Δt*quadfbasis((x)->f(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      end
    elseif alg == "ImplicitEuler"
      methodType = "Implicit"
      K = eye(N) + Δt*Minv*D*A #D okay since numVar = 1 for linear
      lhs = K[freeNode,freeNode]
      if stochastic
        rhs(u,i,dW) = u[freeNode,:] + (Minv*Δt*quadfbasis((x)->f(x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:] +
                    (√Δt.*dW.*Minv*quadfbasis((x)->σ(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      else #Deterministic
        rhs(u,i) = u[freeNode,:] + (Minv*Δt*quadfbasis((x)->f(x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      end
    elseif alg == "CrankNicholson"
      methodType = "Implicit"
      Km = eye(N) - Δt*Minv*D*A/2 #D okay since numVar = 1 for linear
      Kp = eye(N) + Δt*Minv*D*A/2 #D okay since numVar = 1 for linear
      lhs = Kp[freeNode,freeNode]
      if stochastic
        rhs(u,i,dW) = Km[freeNode,freeNode]*u[freeNode,:] + (Minv*Δt*quadfbasis((x)->f(x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:] +
                    (√Δt.*dW.*Minv*quadfbasis((x)->σ(x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      else #Deterministic
        rhs(u,i) = Km[freeNode,freeNode]*u[freeNode,:] + (Minv*Δt*quadfbasis((x)->f(x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      end
    end
  else #Nonlinear Algorithms
    if alg == "Euler"
      methodType = "Explicit"
      if femMesh.μ>=0.5
        warn("Euler method chosen but μ>=.5 => Unstable. Results may be wrong.")
      end
      K = eye(N) - Δt*Minv*A
      if stochastic
        rhs(u,i,dW) = u[freeNode,:] - D.*(Δt*Minv[freeNode,freeNode]*A[freeNode,freeNode]*u[freeNode,:]) + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:] +
                    (√Δt.*dW.*Minv*quadfbasis((u,x)->σ(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      else #Deterministic
        function rhs(u,i)
          u[freeNode,:] - D.*(Δt*Minv[freeNode,freeNode]*A[freeNode,freeNode]*u[freeNode,:]) + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
        end
      end
    elseif alg == "SemiImplicitEuler" #Incorrect for system with different diffusions
      methodType = "Implicit"
      Dinv = D.^(-1)
      K = eye(N) + Δt*Minv*A
      lhs = K[freeNode,freeNode]
      if stochastic
        rhs(u,i,dW) = u[freeNode,:] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:] +
                    (√Δt.*dW.*Minv*quadfbasis((u,x)->σ(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      else #Deterministic
        rhs(u,i) = u[freeNode,:] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      end
    elseif alg == "SemiImplicitCrankNicholson" #Incorrect for system with different diffusions
      methodType = "Implicit"
      Km = eye(N) - Δt*Minv*A/2
      Kp = eye(N) + Δt*Minv*A/2
      lhs = Kp[freeNode,freeNode]
      if stochastic
        rhs(u,i,dW) = Km[freeNode,freeNode]*u[freeNode,:] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:] +
                    (√Δt.*dW.*Minv*quadfbasis((u,x)->σ(u,x,(i-1)*Δt),(x)->gD(x,(i-1)*Δt),(x)->gN(x,(i-1)*Δt),
                                A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      else #Deterministic
        rhs(u,i) = Km[freeNode,freeNode]*u[freeNode,:] + (Minv*Δt*quadfbasis((u,x)->f(u,x,(i-.5)*Δt),(x)->gD(x,(i-.5)*Δt),(x)->gN(x,(i-.5)*Δt),
                    A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
      end
    elseif alg == "ImplicitEuler" # Does this have an issue?
      methodType = "NonlinearSolve"
      if stochastic
        function rhs!(u,resid,dW,uOld,i)
          u = reshape(u,N,numVars)
          resid = reshape(resid,N,numVars)
          resid[freeNode,:] = u[freeNode,:] - uOld[freeNode,:] + D.*(Δt*Minv*A[freeNode,freeNode]*u[freeNode,:]) - (Minv*Δt*quadfbasis((u,x)->f(u,x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:] -(√Δt.*dW.*Minv*quadfbasis((u,x)->σ(u,x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),
                      A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
          u = vec(u)
          resid = vec(resid)
        end
      else #Deterministic
        function rhs!(u,resid,uOld,i)
          u = reshape(u,N,numVars)
          uOld = reshape(uOld,N,numVars)
          resid = reshape(resid,N,numVars)
          resid[freeNode,:] = u[freeNode,:] - uOld[freeNode,:] + D.*(Δt*Minv*A[freeNode,freeNode]*u[freeNode,:]) -
          (Minv*Δt*quadfbasis((u,x)->f(u,x,(i)*Δt),(x)->gD(x,(i)*Δt),(x)->gN(x,(i)*Δt),A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars))[freeNode,:]
          u = vec(u)
          resid = vec(resid)
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
  if stochastic
    rands = getNoise(u,node,elem,noiseType=noiseType)
  end
  if methodType == "NonlinearSolve"
    uOld = similar(vec(u))
  end
  #Heat Equation Loop
  for i=1:femMesh.numIters
    t += Δt
    if methodType == "Implicit"
      if stochastic
        dW = next(rands)
        if solver == "Direct" || solver == "Cholesky" || solver == "QR" || solver == "LU" || solver == "SVD"
          u[freeNode,:] = lhs\rhs(u,i,dW)
        elseif solver == "CG"
          u[freeNode],ch = cg!(u[freeNode],lhs,rhs(u,i,dW)) # Requires Vector, need to change rhs
        elseif solver == "GMRES"
          u[freeNode],ch = gmres!(u[freeNode],lhs,rhs(u,i,dW)) # Requires Vector, need to change rhs
        end
      else #Deterministic
        if solver == "Direct" || solver == "Cholesky" || solver == "QR" || solver == "LU" || solver == "SVD"
          u[freeNode,:] = lhs\rhs(u,i)
        elseif solver == "CG"
          u[freeNode],ch = cg!(u[freeNode],lhs,(u,i)->vec(rhs(u,i))) # Requires Vector, need to change rhs
        elseif solver == "GMRES"
          u[freeNode],ch = gmres!(u[freeNode],lhs,(u,i)->vec(rhs(u,i))) # Requires Vector, need to change rhs
        end
      end
    elseif methodType == "Explicit"
      if stochastic
        dW = next(rands)
        u[freeNode,:] = rhs(u,i,dW)
      else
        u[freeNode,:] = rhs(u,i)
      end
    elseif methodType == "NonlinearSolve"
      u = vec(u)
      uOld = copy(u)
      if stochastic
        dW = next(rands)
        nlres = nlsolve((u,resid)->rhs!(u,resid,dW,uOld,i),uOld,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
      else
        nlres = nlsolve((u,resid)->rhs!(u,resid,uOld,i),uOld,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
      end
      u = nlres.zero
      if numVars > 1
        u = reshape(u,N,numVars)
      end
    end
    u[bdNode] = gD(node,i*Δt)[bdNode]
    if fullSave && i%saveSteps==0
      push!(uFull,u)
      push!(tFull,t)
    end
    atomLoaded ? Main.Atom.progress(i/femMesh.numIters) : nothing #Use Atom's progressbar if loaded
  end
  if knownSol #True Solution exists
    if fullSave
      println("here")
      uFull = copy(uFull)
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
quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear;gNquad𝒪=2)

Performs the order 2 quadrature to calculate the vector from the term ``<f,v>``.
"""
function quadfbasis(f,gD,gN,A,u,node,elem,area,bdNode,mid,N,Dirichlet,Neumann,isLinear,numVars;gNquad𝒪=2)
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
  b = Array{eltype(bt1)}(N,numVars) #size(bt1,2) == numVars
  for i = 1:numVars
    b[:,i] = accumarray(vec(elem),vec([bt1[:,i];bt2[:,i];bt3[:,i]]))
  end

  if(!isempty(Dirichlet))
    uz = zeros(u)
    uz[bdNode,:] = gD(node[bdNode,:])
    b = b-A*uz
  end

  if(!isempty(Neumann))
    el = sqrt(float(sum((node[Neumann[:,1],:] - node[Neumann[:,2],:]).^2,2)))
    λgN,ωgN = quadpts1(gNquad𝒪)
    ϕgN = λgN                # linear bases
    nQuadgN = size(λgN,1)
    ge = zeros(size(Neumann,1),2,numVars)

    for pp = 1:nQuadgN
        # quadrature points in the x-y coordinate
        ppxy = λgN[pp,1]*node[Neumann[:,1],:] +
               λgN[pp,2]*node[Neumann[:,2],:]
        gNp = gN(ppxy)
        for igN = 1:2, var = 1:numVars
            ge[:,igN,var] = ge[:,igN,var] + vec(ωgN[pp]*ϕgN[pp,igN]*gNp[:,var])
        end
    end
    ge = ge.*repeat(el,outer=[1,2,numVars]) # tuple in v0.5?

    for i=1:numVars
      b[:,i] = b[:,i] + accumarray(vec(Neumann), vec(ge[:,i]),[N,1])
    end
  end
  if numVars == 1
    b = vec(b)
  end
  return(b)
end
