immutable FEMHeatIntegrator{T1,T2,T3}
  N::Int
  NT::Int
  Δt::Float64
  t::Number
  Minv::AbstractArray
  D#::AbstractArray
  A::AbstractArray
  freenode::AbstractArray
  f::Function
  gD::Function
  gN::Function
  u::AbstractArray
  node::AbstractArray
  elem::AbstractArray
  area::AbstractArray
  bdnode::AbstractArray
  mid::AbstractArray
  dirichlet::AbstractArray
  neumann::AbstractArray
  islinear::Bool
  numvars::Int
  sqrtΔt::Float64
  σ::Function
  noisetype::Symbol
  numiters::Int
  save_timeseries::Bool
  timeseries#::Vector{uType}
  ts::AbstractArray
  atomloaded::Bool
  solver::Symbol
  autodiff::Bool
  method::Symbol
  show_trace::Bool
  iterations::Int
  timeseries_steps::Int
  progressbar::Bool
  progress_steps::Int
end

@def femheat_footer begin
  u[bdnode] = gD(i*Δt,node)[bdnode]
  if save_timeseries && i%timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
  end
  (atomloaded && progressbar && i%progress_steps==0) ? Main.Atom.progress(i/numiters) : nothing #Use Atom's progressbar if loaded
end

@def femheat_deterministicimplicitlinearsolve begin
  if solver==:Direct || solver==:Cholesky || solver==:QR || solver==:LU || solver==:SVD
    u[freenode,:] = lhs\(Dinv.*rhs(i,u))
  elseif solver==:CG
    for j=1:size(u,2)
      u[freenode,j],ch = cg!(u[freenode,j],lhs,Dinv.*rhs(i,u)[:,j]) # Requires Vector, need to change rhs
    end
  elseif solver==:GMRES
    for j=1:size(u,2)
      u[freenode,j],ch = gmres!(u[freenode,j],lhs,Dinv.*rhs(i,u)[:,j]) # Requires Vector, need to change rhs
    end
  end
end

@def femheat_stochasticimplicitlinearsolve begin
  dW = next(rands)
  if solver==:Direct || solver==:Cholesky || solver==:QR || solver==:LU || solver==:SVD
    u[freenode,:] = lhs\(Dinv.*rhs(i,u,dW))
  elseif solver==:CG
    for j=1:size(u,2)
      u[freenode,j],ch = cg!(u[freenode,j],lhs,Dinv.*rhs(i,u,dW)[:,j]) # Requires Vector, need to change rhs
    end
  elseif solver==:GMRES
    for j=1:size(u,2)
      u[freenode,j],ch = gmres!(u[freenode,j],lhs,Dinv.*rhs(i,u,dW)[:,j]) # Requires Vector, need to change rhs
    end
  end
end

@def femheat_implicitpreamble begin
  Dinv = D.^(-1)
  if solver==:Cholesky
    lhs = cholfact(lhs) # Requires positive definite, may be violated
  elseif solver==:LU
    lhs = lufact(lhs)
  elseif solver==:QR
    lhs = qrfact(lhs) #More stable, slower than LU
  elseif solver==:SVD
    lhs = svdfact(lhs)
  end
end

@def femheat_nonlinearsolvestochasticloop begin
  u = vec(u)
  uOld = copy(u)
  dW = next(rands)
  nlres = NLsolve.nlsolve((u,resid)->rhs!(i,u,resid,dW,uOld),uOld,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
  u = nlres.zero
  if numvars > 1
    u = reshape(u,N,numvars)
  end
end

@def femheat_nonlinearsolvedeterministicloop begin
  u = vec(u)
  uOld = copy(u)
  nlres = NLsolve.nlsolve((u,resid)->rhs!(i,u,resid,uOld),uOld,autodiff=autodiff,method=method,show_trace=show_trace,iterations=iterations)
  u = nlres.zero
  if numvars > 1
    u = reshape(u,N,numvars)
  end
end

@def femheat_nonlinearsolvepreamble begin
  initialize_backend(:NLsolve)
  uOld = similar(vec(u))
end

@def femheat_deterministicpreamble begin
  @unpack N,NT,Δt,t,Minv,D,A,freenode,f,gD,gN,u,node,elem,area,bdnode,mid,dirichlet,neumann,islinear,numvars,numiters,save_timeseries,timeseries,ts,atomloaded,solver,autodiff,method,show_trace,iterations,timeseries_steps,progressbar,progress_steps = integrator
end

@def femheat_stochasticpreamble begin
  @unpack sqrtΔt,σ,noisetype = integrator
  rands = getNoise(u,node,elem,noisetype=noisetype)
end

function femheat_solve(integrator::FEMHeatIntegrator{:linear,:Euler,:deterministic})
  @femheat_deterministicpreamble
  K = eye(N) - Δt*Minv*D*A #D okay since numVar = 1 for linear
  @inbounds for i=1:numiters
    u[freenode,:] = K[freenode,freenode]*u[freenode,:] + (Minv*Δt*quadfbasis((x)->f(t,x),(x)->gD(t,x),(x)->gN(t,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    t += Δt
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:linear,:Euler,:stochastic})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  K = eye(N) - Δt*Minv*D*A #D okay since numVar = 1 for linear
  @inbounds for i=1:numiters
    dW = next(rands)
    u[freenode,:] = K[freenode,freenode]*u[freenode,:] + (Minv*Δt*quadfbasis((x)->f(t,x),(x)->gD(t,x),(x)->gN(t,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
                (sqrtΔt.*dW.*Minv*quadfbasis((x)->σ(t,x),(x)->gD(t,x),(x)->gN(t,x),
                            A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    t += Δt
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:Euler,:deterministic})
  @femheat_deterministicpreamble
  @inbounds for i=1:numiters
    u[freenode,:] = u[freenode,:] - D.*(Δt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) + (Minv*Δt*quadfbasis((x,u)->f(t,x,u),(x)->gD(t,x),(x)->gN(t,x),
            A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    t += Δt
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:Euler,:stochastic})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    u[freenode,:] = u[freenode,:] - D.*(Δt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) + (Minv*Δt*quadfbasis((x,u)->f(t,x,u),(x)->gD(t,x),(x)->gN(t,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
                (sqrtΔt.*dW.*(Minv*quadfbasis((x,u)->σ(t,x,u),(x)->gD(t,x),(x)->gN(t,x),
                            A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
    t += Δt
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:linear,:ImplicitEuler,:stochastic})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  K = eye(N) + Δt*Minv*D*A #D okay since numVar = 1 for linear
  lhs = K[freenode,freenode]
  rhs(i,u,dW) = u[freenode,:] + (Minv*Δt*quadfbasis((x)->f((i)*Δt,x),(x)->gD((i)*Δt,x),(x)->gN((i)*Δt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtΔt.*dW.*(Minv*quadfbasis((x)->σ((i-1)*Δt,x),(x)->gD((i-1)*Δt,x),(x)->gN((i-1)*Δt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += Δt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:linear,:ImplicitEuler,:deterministic})
  @femheat_deterministicpreamble
  K = eye(N) + Δt*Minv*D*A #D okay since numVar = 1 for linear
  lhs = K[freenode,freenode]
  rhs(i,u) = u[freenode,:] + (Minv*Δt*quadfbasis((x)->f((i)*Δt,x),(x)->gD((i)*Δt,x),(x)->gN((i)*Δt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += Δt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:linear,:CrankNicholson,:stochastic})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  Km = eye(N) - Δt*Minv*D*A/2 #D okay since numVar = 1 for linear
  Kp = eye(N) + Δt*Minv*D*A/2 #D okay since numVar = 1 for linear
  lhs = Kp[freenode,freenode]
  rhs(i,u,dW) = Km[freenode,freenode]*u[freenode,:] + (Minv*Δt*quadfbasis((x)->f((i-.5)*Δt,x),(x)->gD((i-.5)*Δt,x),(x)->gN((i-.5)*Δt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtΔt.*dW.*(Minv*quadfbasis((x)->σ((i-1)*Δt,x),(x)->gD((i-1)*Δt,x),(x)->gN((i-1)*Δt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += Δt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:linear,:CrankNicholson,:deterministic})
  @femheat_deterministicpreamble
  Km = eye(N) - Δt*Minv*D*A/2 #D okay since numVar = 1 for linear
  Kp = eye(N) + Δt*Minv*D*A/2 #D okay since numVar = 1 for linear
  lhs = Kp[freenode,freenode]
  rhs(i,u) = Km[freenode,freenode]*u[freenode,:] + (Minv*Δt*quadfbasis((x)->f((i-.5)*Δt,x),(x)->gD((i-.5)*Δt,x),(x)->gN((i-.5)*Δt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += Δt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:SemiImplicitEuler,:deterministic}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  Dinv = D.^(-1)
  K = eye(N) + Δt*Minv*A
  lhs = K[freenode,freenode]
  rhs(i,u) = u[freenode,:] + (Minv*Δt*quadfbasis((x,u)->f((i)*Δt,x,u),(x)->gD((i)*Δt,x),(x)->gN((i)*Δt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += Δt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:SemiImplicitEuler,:stochastic}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  Dinv = D.^(-1)
  K = eye(N) + Δt*Minv*A
  lhs = K[freenode,freenode]
  rhs(i,u,dW) = u[freenode,:] + (Minv*Δt*quadfbasis((x,u)->f((i)*Δt,x,u),(x)->gD((i)*Δt,x),(x)->gN((i)*Δt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtΔt.*dW.*(Minv*quadfbasis((x,u)->σ((i-1)*Δt,x,u),(x)->gD((i-1)*Δt,x),(x)->gN((i-1)*Δt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += Δt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:SemiImplicitCrankNicholson,:deterministic}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  Dinv = D.^(-1)
  Km = eye(N) - Δt*Minv*A/2
  Kp = eye(N) + Δt*Minv*A/2
  lhs = Kp[freenode,freenode]
  rhs(i,u) = Km[freenode,freenode]*u[freenode,:] + (Minv*Δt*quadfbasis((x,u)->f((i-.5)*Δt,x,u),(x)->gD((i-.5)*Δt,x),(x)->gN((i-.5)*Δt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    t += Δt
    @femheat_deterministicimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:SemiImplicitCrankNicholson,:stochastic}) #Incorrect for system with different diffusions
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  Dinv = D.^(-1)
  Km = eye(N) - Δt*Minv*A/2
  Kp = eye(N) + Δt*Minv*A/2
  lhs = Kp[freenode,freenode]
  rhs(i,u,dW) = Km[freenode,freenode]*u[freenode,:] + (Minv*Δt*quadfbasis((x,u)->f((i-.5)*Δt,x,u),(x)->gD((i-.5)*Δt,x),(x)->gN((i-.5)*Δt,x),
              A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] +
              (sqrtΔt.*dW.*(Minv*quadfbasis((x,u)->σ((i-1)*Δt,x,u),(x)->gD((i-1)*Δt,x),(x)->gN((i-1)*Δt,x),
                          A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
  @femheat_implicitpreamble
  @inbounds for i=1:numiters
    dW = next(rands)
    t += Δt
    @femheat_stochasticimplicitlinearsolve
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:ImplicitEuler,:deterministic})
  @femheat_deterministicpreamble
  function rhs!(i,u,resid,uOld)
    u = reshape(u,N,numvars)
    uOld = reshape(uOld,N,numvars)
    resid = reshape(resid,N,numvars)
    resid[freenode,:] = u[freenode,:] - uOld[freenode,:] + D.*(Δt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) -
    (Minv*Δt*quadfbasis((x,u)->f((i)*Δt,x,u),(x)->gD((i)*Δt,x),(x)->gN((i)*Δt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:]
    u = vec(u)
    resid = vec(resid)
  end
  @femheat_nonlinearsolvepreamble
  @inbounds for i=1:numiters
    t += Δt
    @femheat_nonlinearsolvedeterministicloop
    @femheat_footer
  end
  u,timeseries,ts
end

function femheat_solve(integrator::FEMHeatIntegrator{:nonlinear,:ImplicitEuler,:stochastic})
  @femheat_deterministicpreamble
  @femheat_stochasticpreamble
  function rhs!(i,u,resid,dW,uOld)
    u = reshape(u,N,numvars)
    resid = reshape(resid,N,numvars)
    resid[freenode,:] = u[freenode,:] - uOld[freenode,:] + D.*(Δt*Minv[freenode,freenode]*A[freenode,freenode]*u[freenode,:]) - (Minv*Δt*quadfbasis((x,u)->f((i)*Δt,x,u),(x)->gD((i)*Δt,x),(x)->gN((i)*Δt,x),A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars))[freenode,:] -(sqrtΔt.*dW.*(Minv*quadfbasis((x,u)->σ((i)*Δt,x,u),(x)->gD((i)*Δt,x),(x)->gN((i)*Δt,x),
                A,u,node,elem,area,bdnode,mid,N,NT,dirichlet,neumann,islinear,numvars)))[freenode,:]
    u = vec(u)
    resid = vec(resid)
  end
  @femheat_nonlinearsolvepreamble
  @inbounds for i=1:numiters
    t += Δt
    @femheat_nonlinearsolvestochasticloop
    @femheat_footer
  end
  u,timeseries,ts
end
