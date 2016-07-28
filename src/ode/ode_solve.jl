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
function solve(prob::ODEProblem,tspan::AbstractArray=[0,1];kwargs...)
  tspan = vec(tspan)
  if tspan[2]-tspan[1]<0 || length(tspan)>2
    error("tspan must be two numbers and final time must be greater than starting time. Aborting.")
  end

  o = KW(kwargs)
  t = tspan[1]
  T = tspan[2]
  o[:t] = t
  o[:T] = tspan[2]
  @unpack prob: f,u₀,knownSol,sol,numVars,sizeu


  if typeof(u₀)<:Number
    uType = typeof(u₀)
  elseif typeof(u₀) <: AbstractArray
    uType = eltype(u₀)
  else
    error("u₀ must be a number or an array")
  end

  u = u₀

  if :alg ∈ keys(o)
    alg = o[:alg]
  else
    alg = :ExplicitRK
  end

  if alg ∈ DIFFERENTIALEQUATIONSJL_ALGORITHMS

    o2 = copy(DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS)
    for (k,v) in o
      o2[k]=v
    end
    o = o2
    Δt = o[:Δt]
    if Δt==0
      Δt = ode_determine_initΔt(u₀,float(tspan[1]),abstol,reltol,internalNorm,f,order)
    end
    if alg ∉ DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS
      o[:adaptive] = false
    else
      Δt = float(Δt)
    end
    tType=typeof(Δt)

    if o[:Δtmax] == nothing
      o[:Δtmax] = tType((tspan[2]-tspan[1])//2)
    end
    if o[:Δtmin] == nothing
      o[:Δtmin] = tType(1//10^(10))
    end

    T = tType(T)
    t = tType(t)
    uFull = GrowableArray(u₀)
    tFull = Vector{tType}(0)
    push!(tFull,t)
    order = DIFFERENTIALEQUATIONSJL_ORDERS[alg]
    if alg==:ExplicitRK
      @unpack o[:tableau]: A,c,α,αEEst,stages,order
    end
    @materialize maxiters,saveSteps,fullSave,adaptive,progressBar,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,tableau = o
    iter = 0
    if alg==:Euler
      u,t,uFull,tFull = ode_euler(f,u,t,Δt,T,iter,maxiters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
    elseif alg==:Midpoint
      u,t,uFull,tFull = ode_midpoint(f,u,t,Δt,T,iter,maxiters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
    elseif alg==:RK4
      u,t,uFull,tFull = ode_rk4(f,u,t,Δt,T,iter,maxiters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
    elseif alg==:ExplicitRK
      u,t,uFull,tFull = ode_explicitrk(f,u,t,Δt,T,iter,maxiters,uFull,tFull,saveSteps,fullSave,A,c,α,αEEst,stages,order,γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar)
    elseif alg==:ImplicitEuler
      u,t,uFull,tFull = ode_impliciteuler(f,u,t,Δt,T,iter,maxiters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu,progressBar)
    elseif alg==:Trapezoid
      u,t,uFull,tFull = ode_trapezoid(f,u,t,Δt,T,iter,maxiters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu,progressBar)
    elseif alg==:Rosenbrock32
      u,t,uFull,tFull = ode_rosenbrock32(f,u,t,Δt,T,iter,maxiters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar,γ)
    end

  elseif alg ∈ ODEINTERFACE_ALGORITHMS

    if typeof(u) <: Number
      u = [u]
    end
    initialize_backend(:ODEInterface)
    dict = buildOptions(o,ODEINTERFACE_OPTION_LIST,ODEINTERFACE_ALIASES,ODEINTERFACE_ALIASES_REVERSED)
    opts = ODEInterface.OptionsODE([Pair(ODEINTERFACE_STRINGS[k],v) for (k,v) in dict]...) #Convert to the strings
    if alg==:dopri5
      tFull,vecuFull,retcode,stats = ODEInterface.odecall(ODEInterface.dopri5,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opts)
    elseif alg==:dop853
      tFull,vecuFull,retcode,stats = ODEInterface.odecall(ODEInterface.dop853,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opts)
    elseif alg==:odex
      tFull,vecuFull,retcode,stats = ODEInterface.odecall(ODEInterface.odex,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opts)
    elseif alg==:seulex
      tFull,vecuFull,retcode,stats = ODEInterface.odecall(ODEInterface.seulex,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opts)
    elseif alg==:radau
      tFull,vecuFull,retcode,stats = ODEInterface.odecall(ODEInterface.radau,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opts)
    elseif alg==:radau5
      tFull,vecuFull,retcode,stats = ODEInterface.odecall(ODEInterface.radau5,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opts)
    end
    t = tFull[end]
    if typeof(u₀)<:AbstractArray
      uFull = GrowableArray(u₀;initvalue=false)
      for i=1:size(vecuFull,1)
        push!(uFull,reshape(vecuFull[i,:]',sizeu))
      end
    else
      uFull = vecuFull
    end
    u = uFull[end]

  elseif alg ∈ ODEJL_ALGORITHMS
    if typeof(u) <: Number
      u = [u]
    end
    # Needs robustness
    o[:T] = float(o[:T])
    o[:t] = float(o[:t])
    t = o[:t]
    initialize_backend(:ODEJL)
    opts = buildOptions(o,ODEJL_OPTION_LIST,ODEJL_ALIASES,ODEJL_ALIASES_REVERSED)

    ode  = ODE.ExplicitODE(t,u,(t,y,dy)->dy[:]=f(y,t)) #not really inplace
    # adaptive==true ? FoA=:adaptive : FoA=:fixed #Currently limied to only adaptive
    FoA = :adaptive
    if alg==:ode23
      stepper = ODE.RKIntegrator{FoA,:rk23}
    elseif alg==:ode45
      stepper = ODE.RKIntegrator{FoA,:dopri5}
    elseif alg==:ode78
      stepper = ODE.RKIntegrator{FoA,:feh78}
    elseif alg==:ode23s
      stepper = ODE.ModifiedRosenbrockIntegrator
    elseif alg==:ode1
      stepper = ODE.RKIntegratorFixed{:feuler}
    elseif alg==:ode2_midpoint
      stepper = ODE.RKIntegratorFixed{:midpoint}
    elseif alg==:ode2_heun
      stepper = ODE.RKIntegratorFixed{:heun}
    elseif alg==:ode4
      stepper = ODE.RKIntegratorFixed{:rk4}
    elseif alg==:ode45_fe
      stepper = ODE.RKIntegrator{FoA,:rk45}
    end
    out = collect(ODE.solve(ode,stepper;opts...))
    uFull = GrowableArray(u₀)
    tFull = Vector{typeof(out[1][1])}(0)
    push!(tFull,t)
    for (t,u,du) in out
      push!(tFull,t)
      if typeof(u₀) <: AbstractArray
        push!(uFull,u)
      else
        push!(uFull,u[1])
      end
    end
    t = tFull[end]
    u = uFull[end]
  end

  if knownSol
    uTrue = sol(u₀,t)
    if o[:fullSave]
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
    if o[:fullSave]
      uFull = copy(uFull)
      return(ODESolution(u,uFull=uFull,tFull=tFull))
    else
      return(ODESolution(u))
    end
  end
end

function buildOptions(o,optionlist,aliases,aliases_reversed)
  dict1 = Dict{Symbol,Any}([Pair(k,o[k]) for k in (keys(o) ∩ optionlist)])
  dict2 = Dict([Pair(aliases_reversed[k],o[k]) for k in (keys(o) ∩ values(aliases))])
  merge(dict1,dict2)
end

function ode_determine_initΔt(u₀,t,abstol,reltol,internalNorm,f,order)
  d₀ = norm(u₀./(abstol+u₀*reltol),internalNorm)
  f₀ = f(u₀,t)
  d₁ = norm(f₀./(abstol+u₀*reltol),internalNorm)
  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    Δt₀ = 1//10^(6)
  else
    Δt₀ = (d₀/d₁)/100
  end
  u₁ = u₀ + Δt₀*f₀
  f₁ = f(u₁,t+Δt₀)
  d₂ = norm((f₁-f₀)./(abstol+u₀*reltol),internalNorm)/Δt₀
  if max(d₁,d₂)<=1//10^(15)
    Δt₁ = max(1//10^(6),Δt₀*1//10^(3))
  else
    Δt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order+1))
  end
  Δt = min(100Δt₀,Δt₁)
end
