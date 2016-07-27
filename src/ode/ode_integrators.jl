@def ode_loopheader begin
  iter += 1
  if iter > maxIters
    warn("Max Iters Reached. Aborting")
    # u = map((x)->oftype(x,NaN),u)
    break
  end
end

@def ode_savevalues begin
  if fullSave && iter%saveSteps==0
    push!(uFull,u)
    push!(tFull,t)
  end
end

@def ode_loopfooter begin
  if adaptive
    standard = abs(1/(γ*EEst))^(1/order)
    if isinf(standard)
        q = qmax
    else
       q = min(qmax,max(standard,eps()))
    end
    if q > 1
      t = t + Δt
      u = utmp
      @ode_savevalues
    end
    Δtpropose = min(Δtmax,q*Δt)
    Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
  else #Not adaptive
    t = t + Δt
    @ode_savevalues
  end
  (progressBar && atomLoaded && iter%progressSteps==0) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded
end

function ode_euler(f::Function,u,t,Δt,T,iter,maxIters,
                    uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
  while t < T
    @ode_loopheader
    u = u + Δt.*f(u,t)
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_midpoint(f::Function,u::Number,t,Δt,T,iter,
                      maxIters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
  halfΔt = Δt/2
  while t < T
    @ode_loopheader
    u = u + Δt.*f(u+halfΔt.*f(u,t),t+halfΔt)
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_midpoint(f::Function,u::AbstractArray,t,Δt,T,iter,
                      maxIters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
  halfΔt = Δt/2
  utilde = similar(u)
  while t < T
    @ode_loopheader
    utilde[:] = u+halfΔt.*f(u,t)
    u = u + Δt.*f(utilde,t+halfΔt)
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_rk4(f::Function,u::Number,t,Δt,T,iter,maxIters,
                uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
  halfΔt = Δt/2
  while t < T
    @ode_loopheader
    k₁ = f(u,t)
    ttmp = t+halfΔt
    k₂ = f(u+halfΔt*k₁,ttmp)
    k₃ = f(u+halfΔt*k₂,ttmp)
    k₄ = f(u+Δt*k₃,t+Δt)
    u = u + Δt*(k₁ + 2k₂ + 2k₃ + k₄)/6
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_rk4(f::Function,u::AbstractArray,t,Δt,T,
                iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
  halfΔt = Δt/2
  k₁ = similar(u)
  k₂ = similar(u)
  k₃ = similar(u)
  k₄ = similar(u)
  while t < T
    @ode_loopheader
    k₁[:] = f(u,t)
    ttmp = t+halfΔt
    k₂[:] = f(u+halfΔt*k₁,ttmp)
    k₃[:] = f(u+halfΔt*k₂,ttmp)
    k₄[:] = f(u+Δt*k₃,t+Δt)
    u = u + Δt*(k₁ + 2k₂ + 2k₃ + k₄)/6
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_explicitrk(f::Function,u::Number,t,Δt,T,iter,maxIters,uFull,tFull,
                        saveSteps,fullSave,A,c,α,αEEst,stages,order,adaptive,
                        abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar)
  ks = Array{typeof(u)}(stages)
  while t < T
    @ode_loopheader
    for i = 1:stages
      utilde = 0
      for j = 1:i-1
        utilde += A[i,j]*ks[j]
      end
      ks[i] = f(u+Δt*utilde,t+c[i]*Δt)
    end
    utilde = α[1]*ks[1]
    for i = 2:stages
      utilde += α[i]*ks[i]
    end
    if adaptive
      utmp = u + Δt*utilde
      uEEst = αEEst[1]*ks[1]
      for i = 2:stages
        uEEst += αEEst[i]*ks[i]
      end
      EEst = norm((utilde-uEEst)./(abstol+u*reltol),internalNorm)
    else
      u = u + Δt*utilde
    end
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_explicitrk(f::Function,u::AbstractArray,t,Δt,T,iter,
                        maxIters,uFull,tFull,saveSteps,fullSave,
                        A,c,α,αEEst,stages,order,adaptive,
                        abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar)
  ks = Array{eltype(u)}(size(u)...,stages)
  while t < T
    @ode_loopheader
    for i = 1:stages
      utilde = zeros(u)
      for j = 1:i-1
        utilde += A[i,j]*ks[..,j]
      end
      ks[..,i] = f(u+Δt*utilde,t+c[i]*Δt)
    end
    utilde = α[1]*ks[..,1]
    for i = 2:stages
      utilde += α[i]*ks[..,i]
    end
    if adaptive
      utmp = u + Δt*utilde
      uEEst = αEEst[1]*ks[..,1]
      for i = 2:stages
        uEEst += αEEst[i]*ks[..,i]
      end
      EEst = norm((utilde-uEEst)./(abstol+u*reltol),internalNorm)
    else
      u = u + Δt*utilde
    end
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_impliciteuler(f::Function,u,t,Δt,T,iter,maxIters,
                            uFull,tFull,saveSteps,fullSave,adaptive,sizeu,progressBar)
  function rhsIE(u,resid,uOld,t,Δt)
    u = reshape(u,sizeu...)
    resid = reshape(resid,sizeu...)
    resid[:] = u - uOld - Δt*f(u,t+Δt)
    vec(u)
    vec(resid)
  end
  while t < T
    @ode_loopheader
    uOld = copy(u)
    u = vec(u)
    nlres = nlsolve((u,resid)->rhsIE(u,resid,uOld,t,Δt),u)
    u = reshape(nlres.zero,sizeu...)
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_trapezoid(f::Function,u,t,Δt,T,iter,maxIters,
                      uFull,tFull,saveSteps,fullSave,adaptive,sizeu,progressBar)
  function rhsTrap(u,resid,uOld,t,Δt)
    u = reshape(u,sizeu...)
    resid = reshape(resid,sizeu...)
    resid[:] = u - uOld - Δt*(f(u,t+Δt)+f(uOld,t))/2
    u = vec(u)
    resid = vec(resid)
  end
  while t < T
    @ode_loopheader
    uOld = copy(u)
    u = vec(u)
    nlres = nlsolve((u,resid)->rhsTrap(u,resid,uOld,t,Δt),u)
    u = reshape(nlres.zero,sizeu...)
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_rosenbrock32(f::Function,u::AbstractArray,t,Δt,T,iter,
                          maxIters,uFull,tFull,saveSteps,fullSave,adaptive,
                          sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar)
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  k₁ = similar(u)
  k₂ = similar(u)
  k₃ = similar(u)
  function vecf(u,t)
    return(vec(f(reshape(u,sizeu...),t)))
  end
  while t < T
    @ode_loopheader
    # Time derivative
    dT = ForwardDiff.derivative((t)->f(u,t),t)
    J = ForwardDiff.jacobian((u)->vecf(u,t),vec(u))
    W = one(J)-Δt*d*J
    f₀ = f(u,t)
    k₁[:] = reshape(W\vec(f₀ + Δt*d*dT),sizeu...)
    f₁ = f(u+Δt*k₁/2,t+Δt/2)
    k₂[:] = reshape(W\vec(f₁-k₁),sizeu...) + k₁
    if adaptive
      utmp = u + Δt*k₂
      f₂ = f(utmp,t+Δt)
      k₃[:] = reshape(W\vec(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T),sizeu...)
      EEst = norm((Δt(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalNorm)
    else
      u = u + Δt*k₂
    end
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end

function ode_rosenbrock32(f::Function,u::Number,t,Δt,T,iter,
                          maxIters,uFull,tFull,saveSteps,fullSave,adaptive,
                          sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar)
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  function vecf(u,t)
    return(vec(f(reshape(u,sizeu...),t)))
  end
  while t < T
    @ode_loopheader
    # Time derivative
    dT = ForwardDiff.derivative((t)->f(u,t),t)
    J = ForwardDiff.jacobian((u)->vecf(u,t),vec(u))
    W = one(J)-Δt*d*J
    f₀ = f(u,t)
    k₁ = reshape(W\vec(f₀ + Δt*d*dT),sizeu...)
    f₁ = f(u+Δt*k₁/2,t+Δt/2)
    k₂ = reshape(W\vec(f₁-k₁),sizeu...) + k₁
    if adaptive
      utmp = u + Δt*k₂
      f₂ = f(utmp,t+Δt)
      k₃ = reshape(W\vec(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T),sizeu...)
      EEst = norm((Δt(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalNorm)
    else
      u = u + Δt*k₂
    end
    @ode_loopfooter
  end
  return u,t,uFull,tFull
end
