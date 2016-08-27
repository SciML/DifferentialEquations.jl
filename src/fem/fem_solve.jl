"""
## Finite Element Poisson Equation Solver

`solve(fem_mesh::FEMmesh,pdeProb::PoissonProblem)`

Takes in a definition for the heat equation ``-Δu = f`` on `fem_mesh` with
functions as defined in `pdeProb`. If `σ` is specified in `pdeProb`, then this
solves the stochastic Poisson equation ``-Δu = f + σdW``.

### Keyword Arguments

* `solver` = Linear solver algorithm. This is the algorithm which is chosen for solving
  the implicit equation `Ax=b`. The default is `LU`. The choices are:

    - `:Direct` = Solves `Ax=b` using `\\`
    - `:CG` = Conjugate-Gradient. Best when the space is very large and ``I ± ΔtM⁻¹A`` is positive definite.
    - `:GMRES` = GMRES. Best when the space is very large and ``I ± ΔtM⁻¹A`` is not positive definite.

* `timeseries_steps` = If `save_timeseries=true`, then this is the number of steps between the saves.
* `autodiff` = Whether or not autodifferentiation (as provided by AutoDiff.jl) is used
  for the nonlinear solving. By default autodiff is false.
* `method` = Method the nonlinear solver uses. Defaults to `:trust_region`.
* `show_trace` = Whether to show the output of the nonlinear solver. Defaults to false.
* `iterations` = Maximum numer of iterations in the nonlinear solver. Defaults to 1000.
"""
function solve(fem_mesh::FEMmesh,prob::PoissonProblem;solver::Symbol=:Direct,autodiff::Bool=false,method=:trust_region,show_trace=false,iterations=1000)
  #Assemble Matrices
  A,M,area = assemblematrix(fem_mesh,lumpflag=true)

  #Unroll some important constants
  @unpack fem_mesh: Δt,bdnode,node,elem,N,NT,freenode,dirichlet,neumann
  @unpack prob: f,Du,f,gD,gN,analytic,knownanalytic,islinear,u₀,numvars,σ,stochastic,noisetype,D

  #Setup f quadrature
  mid = Array{Float64}(size(node[vec(elem[:,2]),:])...,3)
  mid[:,:,1] = (node[vec(elem[:,2]),:]+node[vec(elem[:,3]),:])/2
  mid[:,:,2] = (node[vec(elem[:,3]),:]+node[vec(elem[:,1]),:])/2
  mid[:,:,3] = (node[vec(elem[:,1]),:]+node[vec(elem[:,2]),:])/2

  #Setup u
  u = u₀(node)
  if numvars==0
    numvars = size(u,2)
    prob.numvars = numvars #Mutate problem to be correct.
    if gD == nothing
      gD=(x)->zeros(size(x,1),numvars)
    end
    prob.gD = gD
    if gN == nothing
      gN=(x)->zeros(size(x,1),numvars)
    end
    prob.gN = gN
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
    prob.D = D
  end
  u[bdnode] = gD(node[bdnode,:])

  #Stochastic Part
  if stochastic
    rands = getNoise(u,node,elem,noisetype=noisetype)
    dW = next(rands)
    rhs = (u) -> quadfbasis(f,gD,gN,A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars) + quadfbasis(σ,gD,gN,A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars).*dW
  else #Not Stochastic
    rhs = (u) -> quadfbasis(f,gD,gN,A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)
  end
  Dinv = D.^(-1)
  #Solve
  if islinear
    if solver==:Direct
      u[freenode,:]=D.*(A[freenode,freenode]\rhs(u)[freenode])
    elseif solver==:CG
      for i = 1:size(u,2)
        u[freenode,i],ch=cg!(u[freenode,i],A[freenode,freenode],Dinv.*rhs(u)[freenode,i]) # Needs diffusion constant
      end
    elseif solver==:GMRES
      for i = 1:size(u,2)
        u[freenode,i],ch=gmres!(u[freenode,i],A[freenode,freenode],Dinv.*rhs(u)[freenode,i]) # Needs diffusion constants
      end
    end
    #Adjust result
    if isempty(dirichlet) #isPureneumann
      patchArea = accumarray(vec(elem),[area;area;area]/3, [N 1])
      uc = sum(u.*patchArea)/sum(area)
      u = u - uc   # Impose integral of u = 0
    end
  else #Nonlinear
    rhs! = (u,resid) -> begin
      u = reshape(u,N,numvars)
      resid = reshape(resid,N,numvars)
      resid[freenode,:]=D.*(A[freenode,freenode]*u[freenode,:])-rhs(u)[freenode,:]
      u = vec(u)
      resid = vec(resid)
    end
    initialize_backend(:NLsolve)
    u = vec(u)
    nlres = NLsolve.nlsolve(rhs!,u,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
    u = nlres.zero
    if numvars > 1
      u = reshape(u,N,numvars)
    end
  end

  #Return
  if knownanalytic # True solution exists
    return(FEMSolution(fem_mesh,u,analytic(node),analytic,Du,prob))
  else #No true solution
    return(FEMSolution(fem_mesh,u,prob))
  end
end

## Evolution Equation Solvers
#Note
#rhs(u,i) = Dm[freenode,freenode]*u[freenode,:] + Δt*f(node,(i-.5)*Δt)[freenode] #Nodel interpolation 1st 𝒪
"""
## Finite Element Heat Equation Solver

`solve(fem_mesh::FEMmesh,pdeProb::HeatProblem)`

Takes in a definition for the heat equation ``u_t = Δu + f`` on `fem_mesh` with
functions as defined in `pdeProb`. If `σ` is specified in `pdeProb`, then this
solves the stochastic heat equation ``u_t = Δu + f + σdW_t``.

### Keyword Arguments

* `alg` = Solution algorithm. Default is :Euler. The choices are:

    - Linear

        * `:Euler` (Explicit)
        * `:ImplicitEuler` (Implicit)
        * `:CrankNicholson` (Implicit)

    - Nonlinear

        * `:Euler` (Explicit)
        * `:ImplicitEuler` (Nonlinear Solve)
        * `:CrankNicholson` (Nonlinear Solve)
        * `:SemiImplicitEuler` (Implicit)
        * `:SemiImplicitCrankNicholson` (Implicit)

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

    - `:Direct` = Solves using `\\` (no factorization). Not recommended.
    - `:Cholesky` = Cholsky decomposition. Only stable of ``I ± ΔtM⁻¹A`` is positive definite.
      This means that this works best when Δt is small. When applicable, this is the fastest.
    - `:LU` = LU-Decomposition. A good mix between fast and stable.
    - `:QR` = QR-Decomposition. Less numerical roundoff error than `LU`, but slightly slower.
    - `:SVD` = SVD-Decomposition. By far the slowest, but the most robust to roundoff error.
    - `:CG` = Conjugate-Gradient. Best when the space is very large and ``I ± ΔtM⁻¹A`` is positive definite.
    - `:GMRES` = GMRES. Best when the space is very large and ``I ± ΔtM⁻¹A`` is not positive definite.

* `save_timeseries` = Makes the algorithm save the output at every `timeseries_steps` timesteps.
  By default save_timeseries is false.
* `timeseries_steps` = If `save_timeseries=true`, then this is the number of steps between the saves.
* `autodiff` = Whether or not autodifferentiation (as provided by AutoDiff.jl) is used
  for the nonlinear solving. By default autodiff is false.
* `method` = Method the nonlinear solver uses. Defaults to `:trust_region`.
* `show_trace` = Whether to show the output of the nonlinear solver. Defaults to false.
* `iterations` = Maximum numer of iterations in the nonlinear solver. Defaults to 1000.
* `progress_steps` = The number of steps between updates of the progress bar. Defaults to 1000.
* `progressbar` = Turns on/off use of the Juno progress bar. Defaults to true. Requires Juno.
"""
function solve(fem_mesh::FEMmesh,prob::HeatProblem;alg::Symbol=:Euler,
  solver::Symbol=:LU,save_timeseries::Bool = false,timeseries_steps::Int = 100,
  autodiff::Bool=false,method=:trust_region,show_trace=false,iterations=1000,
  progress_steps::Int=1000,progressbar::Bool=true)
  #Assemble Matrices
  A,M,area = assemblematrix(fem_mesh,lumpflag=true)

  #Unroll some important constants
  @unpack fem_mesh: Δt,bdnode,node,elem,N,NT,freenode,dirichlet,neumann
  @unpack prob: f,u₀,Du,gD,gN,analytic,knownanalytic,islinear,numvars,σ,stochastic,noisetype,D

  #Note if Atom is loaded for progress
  atomloaded = isdefined(Main,:Atom)

  #Set Initial
  u = copy(u₀(node))
  if numvars==0
    numvars = size(u,2)
    prob.numvars = numvars #Mutate problem to be correct.
    if gD == nothing
      gD=(t,x)->zeros(size(x,1),numvars)
    end
    prob.gD = gD
    if gN == nothing
      gN=(t,x)->zeros(size(x,1),numvars)
    end
    prob.gN = gN
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
    prob.D = D
  end
  t = 0
  sqrtΔt= sqrt(Δt)
  #Setup f quadraturef
  mid = Array{Float64}(size(node[vec(elem[:,2]),:])...,3)
  mid[:,:,1] = (node[vec(elem[:,2]),:]+node[vec(elem[:,3]),:])/2
  mid[:,:,2] = (node[vec(elem[:,3]),:]+node[vec(elem[:,1]),:])/2
  mid[:,:,3] = (node[vec(elem[:,1]),:]+node[vec(elem[:,2]),:])/2

  islinear ? linearity=:linear : linearity=:nonlinear
  stochastic ? stochasticity=:stochastic : stochasticity=:deterministic
  #Setup timeseries

  timeseries = Vector{typeof(u)}(0)
  push!(timeseries,u)
  ts = Float64[t]

  if alg==:Euler && fem_mesh.μ>=0.5
    warn("Euler method chosen but μ>=.5 => Unstable. Results may be wrong.")
  end

  #Setup for Calculations
  Minv = sparse(inv(M)) #sparse(Minv) needed until update

  #Heat Equation Loop
  u,timeseres,ts=femheat_solve(FEMHeatIntegrator{linearity,alg,stochasticity}(N,NT,Δt,t,Minv,D,A,freenode,f,gD,gN,u,node,elem,area,bdnode,mid,dirichlet,neumann,islinear,numvars,sqrtΔt,σ,noisetype,fem_mesh.numiters,save_timeseries,timeseries,ts,atomloaded,solver,autodiff,method,show_trace,iterations,timeseries_steps,progressbar,progress_steps))

  (atomloaded && progressbar ) ? Main.Atom.progress(1) : nothing #Use Atom's progressbar if loaded

  if knownanalytic #True Solution exists
    if save_timeseries
      timeseries = FEMSolutionTS(timeseries,numvars)
      return(FEMSolution(fem_mesh,u,analytic(fem_mesh.T,node),analytic,Du,timeseries,ts,prob))
    else
      return(FEMSolution(fem_mesh,u,analytic(fem_mesh.T,node),analytic,Du,prob))
    end
  else #No true solution
    if save_timeseries
      timeseries = FEMSolutionTS(timeseries,numvars)
      return(FEMSolution(fem_mesh,u,timeseries,ts,prob))
    else
      return(FEMSolution(fem_mesh,u,prob))
    end
  end
end

"""
quadfbasis(f,gD,gN,A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars;gNquad𝒪=2)

Performs the order 2 quadrature to calculate the vector from the term ``<f,v>`` for linear elements.
"""
function quadfbasis(f,gD,gN,A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars;gNquad𝒪=2)
  b = zeros(u) #size(bt1,2) == numvars
  if islinear
    bt1 = area.*(f(mid[:,:,2])+f(mid[:,:,3]))/6
    bt2 = area.*(f(mid[:,:,3])+f(mid[:,:,1]))/6
    bt3 = area.*(f(mid[:,:,1])+f(mid[:,:,2]))/6
  else
    u1 = (u[vec(elem[:,2]),:]+u[vec(elem[:,3]),:])/2
    u2 = (u[vec(elem[:,3]),:]+u[vec(elem[:,1]),:])/2
    u3 = (u[vec(elem[:,1]),:]+u[vec(elem[:,2]),:])/2
    bt1 = area.*(f(mid[:,:,2],u2)+f(mid[:,:,3],u3))/6
    bt2 = area.*(f(mid[:,:,3],u3)+f(mid[:,:,1],u1))/6
    bt3 = area.*(f(mid[:,:,1],u1)+f(mid[:,:,2],u2))/6
  end

  for i = 1:numvars # accumarray the bt's
    for j = 1:NT
      b[elem[j,1],i] += bt1[j,i]
    end
    for j = 1:NT
      b[elem[j,2],i] += bt2[j,i]
    end
    for j = 1:NT
      b[elem[j,3],i] += bt3[j,i]
    end
  end

  if(!isempty(dirichlet))
    uz = zeros(u)
    uz[bdnode,:] = gD(node[bdnode,:])
    b = b-A*uz
  end

  if(!isempty(neumann))
    el = sqrt(float(sum((node[neumann[:,1],:] - node[neumann[:,2],:]).^2,2)))
    λgN,ωgN = quadpts1(gNquad𝒪)
    ϕgN = λgN                # linear bases
    nQuadgN = size(λgN,1)
    ge = zeros(size(neumann,1),2,numvars)

    for pp = 1:nQuadgN
        # quadrature points in the x-y coordinate
        ppxy = λgN[pp,1]*node[neumann[:,1],:] +
               λgN[pp,2]*node[neumann[:,2],:]
        gNp = gN(ppxy)
        for igN = 1:2, var = 1:numvars
            ge[:,igN,var] = ge[:,igN,var] + vec(ωgN[pp]*ϕgN[pp,igN]*gNp[:,var])
        end
    end
    ge = ge.*repeat(el,outer=[1,2,numvars]) # tuple in v0.5?

    for i=1:numvars
      b[:,i] = b[:,i] + accumarray(vec(neumann), vec(ge[:,i]),[N,1])
    end
  end
  if numvars == 1
    b = vec(b)
  end
  return(b)
end
