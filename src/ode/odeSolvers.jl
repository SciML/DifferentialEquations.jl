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

end

function ode_euler(f::Function,u,t,Δt,T,iter,maxIters,
                    uFull,tFull,saveSteps,fullSave,adaptive)
  while t < T
    @ode_loopheader
    u = u + Δt.*f(u,t)
    @ode_loopfooter
  end
  return f,u,t,Δt,uFull,tFull
end

function ode_midpoint(f::Function,u::Number,t,Δt,T,iter,
                      maxIters,uFull,tFull,saveSteps,fullSave,adaptive)
  halfΔt = Δt/2
  while t < T
    @ode_loopheader
    u = u + Δt.*f(u+halfΔt.*f(u,t),t+halfΔt)
    @ode_loopfooter
  end
  return f,u,t,Δt,uFull,tFull
end

function ode_midpoint(f::Function,u::AbstractArray,t,Δt,T,iter,
                      maxIters,uFull,tFull,saveSteps,fullSave,adaptive)
  halfΔt = Δt/2
  utilde = similar(u)
  while t < T
    @ode_loopheader
    utilde[:] = u+halfΔt.*f(u,t)
    u = u + Δt.*f(utilde,t+halfΔt)
    @ode_loopfooter
  end
  return f,u,t,Δt,uFull,tFull
end

function ode_rk4(f::Function,u::Number,t,Δt,T,iter,maxIters,
                uFull,tFull,saveSteps,fullSave,adaptive)
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
  return f,u,t,Δt,uFull,tFull
end

function ode_rk4(f::Function,u::AbstractArray,t,Δt,T,
                iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive)
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
  return f,u,t,Δt,uFull,tFull
end

function ode_explicitrk(f::Function,u::Number,t,Δt,T,iter,maxIters,uFull,tFull,
                        saveSteps,fullSave,A,c,α,αEEst,stages,order,adaptive,
                        abstol,reltol,qmax,Δtmax,Δtmin,internalNorm)
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
  return f,u,t,Δt,uFull,tFull
end

function ode_explicitrk(f::Function,u::AbstractArray,t,Δt,T,iter,
                        maxIters,uFull,tFull,saveSteps,fullSave,
                        A,c,α,αEEst,stages,order,adaptive,
                        abstol,reltol,qmax,Δtmax,Δtmin,internalNorm)
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
  return f,u,t,Δt,uFull,tFull
end

function ode_impliciteuler(f::Function,u,t,Δt,T,iter,maxIters,
                            uFull,tFull,saveSteps,fullSave,adaptive,sizeu)
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
  return f,u,t,Δt,uFull,tFull
end

function ode_trapezoid(f::Function,u,t,Δt,T,iter,maxIters,
                      uFull,tFull,saveSteps,fullSave,adaptive,sizeu)
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
  return f,u,t,Δt,uFull,tFull
end

function ode_rosenbrock32(f::Function,u::AbstractArray,t,Δt,T,iter,
                          maxIters,uFull,tFull,saveSteps,fullSave,adaptive,
                          sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm)
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
  return f,u,t,Δt,uFull,tFull
end

function ode_rosenbrock32(f::Function,u::Number,t,Δt,T,iter,
                          maxIters,uFull,tFull,saveSteps,fullSave,adaptive,
                          sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm)
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
  return f,u,t,Δt,uFull,tFull
end

"""
solve(prob::ODEProblem,Δt,T)

Solves the ODE defined by prob with initial Δt on the time interval [0,T].

### Keyword Arguments

* fullSave: Saves the result at every saveSteps steps. Default is false.
saveSteps: If fullSave is true, then the output is saved every saveSteps steps.
* alg: String which defines the solver algorithm. Defult is "RK4". Possibilities are:
  * "Euler" - The canonical forward Euler method.
  * "Midpoint" - The second order midpoint method.
  * "RK4" - The canonical Runge-Kutta Order 4 method.
  * "ExplicitRK" - A general Runge-Kutta solver which takes in a tableau. Can be adaptive.
  * "ImplicitEuler" - A 1st order implicit solver. Unconditionally stable.
  * "Trapezoid" - A second order unconditionally stable implicit solver. Good for highly stiff.
  * "Rosenbrock32" - A fast solver which is good for stiff equations.
* tableau - Takes in an object which defines a tableau. Default is Dormand-Prince 4/5.
* adaptive - Turns on adaptive timestepping for appropriate methods. Default is false.
* tol - The error tolerance of the adaptive method. Default is 1e-4.
* γ - The risk-factor γ in the q equation for adaptive timestepping. Default is 2.
* qmax - Defines the maximum value possible for the adaptive q. Default is 10.
"""
function solve(prob::ODEProblem,tspan::AbstractArray=[0,1];Δt::Number=0,
              fullSave::Bool = false,saveSteps::Int = 1,alg::Symbol=:RK4,
              tableau=DEFAULT_TABLEAU,adaptive=false,γ=2.0,
              abstol=nothing,reltol=nothing,qmax=4,maxIters::Int = round(Int,1e9),
              Δtmax=nothing,Δtmin=nothing,tType=typeof(Δt),internalNorm = 2)
  tspan = vec(tspan)
  if tspan[2]-tspan[1]<0 || length(tspan)>2
    error("tspan must be two numbers and final time must be greater than starting time. Aborting.")
  end

  @unpack prob: f,u₀,knownSol,sol, numVars, sizeu
  if typeof(u₀)<:Number
    uType = typeof(u₀)
  elseif typeof(u₀) <: AbstractArray
    uType = eltype(u₀)
  else
    error("u₀ must be a number or an array")
  end

  if Δtmax == nothing
    Δtmax = tType((tspan[2]-tspan[1])//2)
  end
  if Δtmin == nothing
    Δtmin = tType(1//10^(10))
  end
  if abstol == nothing
    abstol = uType(1//10^8)
  end
  if reltol == nothing
    reltol = uType(1//10^6)
  end

  T = tType(tspan[2])
  t = tType(tspan[1])
  u = u₀

  @unpack tableau: A,c,α,αEEst,stages,order

  uFull = GrowableArray(u)
  tFull = Vector{tType}(0)
  push!(tFull,t)

  iter = 0

  if Δt == 0
    d₀ = norm(u₀./(abstol+u*reltol),internalNorm)
    f₀ = f(u₀,t)
    d₁ = norm(f₀./(abstol+u*reltol),internalNorm)
    if d₀ < 1//10^(5) || d₁ < 1//10^(5)
      Δt₀ = 1//10^(6)
    else
      Δt₀ = (d₀/d₁)/100
    end
    u₁ = u₀ + Δt₀*f₀
    f₁ = f(u₁,t+Δt₀)
    d₂ = norm((f₁-f₀)./(abstol+u*reltol),internalNorm)/Δt₀
    if max(d₁,d₂)<=1//10^(15)
      Δt₁ = max(1//10^(6),Δt₀*1//10^(3))
    else
      if !isdefined(Main,:order)
        order = 1 #Convervative choice
      end
      Δt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order+1))
    end
    Δt = min(100Δt₀,Δt₁)
  end

  if alg==:Euler
    f,u,t,Δt,uFull,tFull = ode_euler(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive)
  elseif alg==:Midpoint
    f,u,t,Δt,uFull,tFull = ode_midpoint(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive)
  elseif alg==:RK4
    f,u,t,Δt,uFull,tFull = ode_rk4(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive)
  elseif alg==:ExplicitRK
    f,u,t,Δt,uFull,tFull = ode_explicitrk(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,A,c,α,αEEst,stages,order,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm)
  elseif alg==:ImplicitEuler
    f,u,t,Δt,uFull,tFull = ode_impliciteuler(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu)
  elseif alg==:Trapezoid
    f,u,t,Δt,uFull,tFull = ode_trapezoid(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu)
  elseif alg==:Rosenbrock32
    f,u,t,Δt,uFull,tFull = ode_rosenbrock32(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm)
  end

  if knownSol
    uTrue = sol(u₀,t)
    if fullSave
      solFull = GrowableArray(sol(u₀,tFull[1]))
      for i in 2:size(uFull,1)
        push!(solFull,sol(u₀,tFull[i]))
      end
      uFull = copy(uFull)
      solFull = copy(solFull)
      return(ODESolution(u,uTrue,uFull=uFull,tFull=tFull,solFull=solFull))
    else
      return(ODESolution(u,uTrue))
    end
  else #No known sol
    if fullSave
      uFull = copy(uFull)
      return(ODESolution(u,uFull=uFull,tFull=tFull))
    else
      return(ODESolution(u))
    end
  end
end
