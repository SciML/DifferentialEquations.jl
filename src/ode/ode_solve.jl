"""
`solve(prob::ODEProblem,tspan)`

Solves the ODE defined by prob on the interval tspan. If not given, tspan defaults to [0,1].

### Keyword Arguments

* `Δt`: Sets the initial stepsize. Defaults to an automatic choice.
* `save_timeseries`: Saves the result at every timeseries_steps steps. Default is true.
* `timeseries_steps`: Denotes how many steps between saving a value for the timeseries. Defaults to 1.
* `tableau`: The tableau for an `:ExplicitRK` algorithm. Defaults to a Dormand-Prince 4/5 method.
* `adaptive` - Turns on adaptive timestepping for appropriate methods. Default is true.
* `γ` - The risk-factor γ in the q equation for adaptive timestepping. Default is .9.
* `timechoicealg` - Chooses the method which is used for making the adaptive timestep choices.
  Default is `:Lund` for Lund stabilization (PI stepsize control). The other
  option is `:Simple` for the standard simple error-based rejection
* `β` - The Lund stabilization β parameter. Defaults are algorithm-dependent.
* `qmax` - Defines the maximum value possible for the adaptive q. Default is 10.
* `ablstol` - Absolute tolerance in adaptive timestepping. Defaults to 1e-3.
* `reltol` - Relative tolerance in adaptive timestepping. Defaults to 1e-6.
* `maxiters` - Maximum number of iterations before stopping. Defaults to 1e9.
* `Δtmax` - Maximum Δt for adaptive timestepping. Defaults to half the timespan.
* `Δtmin` - Minimum Δt for adaptive timestepping. Defaults to 1e-10.
* `autodiff` - Turns on/off the use of autodifferentiation (via ForwardDiff) in the
  implicit solvers which use `NLsolve`. Default is true.
* `internalnorm` - The norm for which error estimates are calculated. Default is 2.
* `progressbar` - Turns on/off the Juno progressbar. Defualt is false.
* `progress_steps` - Numbers of steps between updates of the progress bar. Default is 1000.

* `alg`: The solver algorithm. Defult is `:DP5`. Note that any keyword
  argument available in the external solvers are accessible via keyword arguemnts. For example,
  for the ODEInterface.jl algorithms, one can specify `SSBETA=0.03` as a keyword argument and it will
  do as it states in the ODEInterface.jl documentation. Common options such as `MAXSS` (max stepsize)
  are aliased to one can use the DifferentialEquations.jl syntax `Δtmax` or `MAXSS`. The possibilities for the solvers are:

For a full list of algorithms, please see the solver documentation.
"""
function solve{uType<:Union{AbstractArray,Number},uEltype<:Number}(prob::ODEProblem{uType,uEltype},tspan::AbstractArray=[0,1];kwargs...)
  tspan = vec(tspan)
  if tspan[end]-tspan[1]<0
    error("final time must be greater than starting time. Aborting.")
  end
  atomloaded = isdefined(Main,:Atom)
  o = KW(kwargs)
  o[:t] = tspan[1]
  o[:Ts] = tspan[2:end]
  @unpack prob: u₀,knownanalytic,analytic,numvars,isinplace


  command_opts = merge(o,DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS)
  # Get the control variables
  @materialize save_timeseries, progressbar = command_opts

  u = copy(u₀)

  if :alg ∈ keys(o)
    alg = o[:alg]
  else
    alg = :DP5 # Default algorithm
  end

  if alg ∈ DIFFERENTIALEQUATIONSJL_ALGORITHMS
    o2 = copy(DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS)
    for (k,v) in o
      o2[k]=v
    end
    o = o2
    Δt = o[:Δt]
    order = DIFFERENTIALEQUATIONSJL_ORDERS[alg]
    adaptiveorder = 0
    if alg ∈ DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS
      adaptiveorder = DIFFERENTIALEQUATIONSJL_ADAPTIVEORDERS[alg]
    end
    if alg==:ExplicitRK || alg==:ExplicitRKVectorized
      @unpack o[:tableau]: order,adaptiveorder
    end
    if !isinplace && typeof(u)<:AbstractArray
      f = (du,u,t) -> (du[:] = prob.f(u,t))
    else
      f = prob.f
    end
    if Δt==0
      Δt = ode_determine_initΔt(u₀,float(tspan[1]),o[:abstol],o[:reltol],o[:internalnorm],f,order)
    end
    if alg ∉ DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS
      o[:adaptive] = false
    else
      if o[:adaptive] == true
        Δt = float(Δt)
      end
    end
    if alg ∈ DIFFERENTIALEQUATIONSJL_IMPLICITALGS
      initialize_backend(:NLsolve)
    end
    tType=typeof(Δt)

    if o[:Δtmax] == nothing
      o[:Δtmax] = tType((tspan[end]-tspan[1]))
    end
    if o[:Δtmin] == nothing
      o[:Δtmin] = tType(1//10^(10))
    end
    if o[:fullnormalize] == true
      normfactor = uEltype(1/length(u))
    else
      normfactor = 1
    end
    if o[:β] == nothing # Use default β
      if alg == :DP5 || alg == :DP5Vectorized
        β = 0.10 # More than Hairer's suggestion
      elseif alg == :DP8 || alg == :DP8Vectorized
        β = 0.07 # More than Hairer's suggestion
      else
        β = 0.4 / order
      end
    else
      β = o[:β]
    end
    fsal = false
    if alg ∈ DIFFERENTIALEQUATIONSJL_FASLALGS
      fsal = true
    elseif alg == :ExplicitRK
      @unpack o[:tableau]: fsal
    end

    Ts = map(tType,o[:Ts])
    t = tType(o[:t])
    timeseries = GrowableArray(u₀)
    ts = Vector{tType}(0)
    push!(ts,t)
    @materialize maxiters,timeseries_steps,save_timeseries,adaptive,progress_steps,abstol,reltol,γ,qmax,qmin,Δtmax,Δtmin,internalnorm,tableau,autodiff, timechoicealg,qoldinit= o
    #@code_warntype  ode_solve(ODEIntegrator{alg,uType,uEltype,ndims(u)+1,tType}(f,u,t,Δt,T,maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,Δtmax,Δtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β,timechoicealg,qoldinit,normfactor,fsal))
    u,t,timeseries,ts = ode_solve(ODEIntegrator{alg,uType,uEltype,ndims(u)+1,tType}(f,u,t,Δt,Ts,maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,Δtmax,Δtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β,timechoicealg,qoldinit,normfactor,fsal))

  elseif alg ∈ ODEINTERFACE_ALGORITHMS
    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    o[:T] = float(o[:T])
    o[:t] = float(o[:t])
    t = o[:t]; T = o[:T]
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = vec(prob.f(reshape(u,sizeu),t)); nothing)
    else
      f! = (t,u,du) -> (prob.f(reshape(du,sizeu),reshape(u,sizeu),t); u = vec(u); du=vec(du); nothing)
    end
    initialize_backend(:ODEInterface)
    o[:RHS_CALLMODE] = ODEInterface.RHS_CALL_INSITU
    dict = buildOptions(o,ODEINTERFACE_OPTION_LIST,ODEINTERFACE_ALIASES,ODEINTERFACE_ALIASES_REVERSED)
    opts = ODEInterface.OptionsODE([Pair(ODEINTERFACE_STRINGS[k],v) for (k,v) in dict]...) #Convert to the strings
    du = similar(u)
    if alg==:dopri5
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dopri5,f!,Float64[t,T],vec(u),opts)
    elseif alg==:dop853
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dop853,f!,Float64[t,T],vec(u),opts)
    elseif alg==:odex
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.odex,f!,Float64[t,T],vec(u),opts)
    elseif alg==:seulex
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.seulex,f!,Float64[t,T],vec(u),opts)
    elseif alg==:radau
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau,f!,Float64[t,T],vec(u),opts)
    elseif alg==:radau5
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau5,f!,Float64[t,T],vec(u),opts)
    end
    t = ts[end]
    if typeof(u₀)<:AbstractArray
      timeseries = GrowableArray(u₀;initvalue=false)
      for i=1:size(vectimeseries,1)
        push!(timeseries,reshape(vectimeseries[i,:]',sizeu))
      end
    else
      timeseries = vectimeseries
    end
    u = timeseries[end]

  elseif alg ∈ ODEJL_ALGORITHMS
    if typeof(u) <: Number
      u = [u]
    end
    # Needs robustness
    o[:T] = float(o[:T])
    o[:t] = float(o[:t])
    t = o[:t]; T = o[:T]
    initialize_backend(:ODEJL)
    opts = buildOptions(o,ODEJL_OPTION_LIST,ODEJL_ALIASES,ODEJL_ALIASES_REVERSED)
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = prob.f(u,t))
    else
      f! = (t,u,du) -> (prob.f(du,u,t))
    end
    ode  = ODE.ExplicitODE(t,u,f!)
    # adaptive==true ? FoA=:adaptive : FoA=:fixed #Currently limied to only adaptive
    FoA = :adaptive
    if alg==:ode23
      solver = ODE.RKIntegrator{FoA,:rk23}
    elseif alg==:ode45
      solver = ODE.RKIntegrator{FoA,:dopri5}
    elseif alg==:ode78
      solver = ODE.RKIntegrator{FoA,:feh78}
    elseif alg==:ode23s
      solver = ODE.ModifiedRosenbrockIntegrator
    elseif alg==:ode1
      solver = ODE.RKIntegratorFixed{:feuler}
    elseif alg==:ode2_midpoint
      solver = ODE.RKIntegratorFixed{:midpoint}
    elseif alg==:ode2_heun
      solver = ODE.RKIntegratorFixed{:heun}
    elseif alg==:ode4
      solver = ODE.RKIntegratorFixed{:rk4}
    elseif alg==:ode45_fe
      solver = ODE.RKIntegrator{FoA,:rk45}
    end
    out = collect(ODE.solve(ode;solver=solver,opts...))
    timeseries = GrowableArray(u₀)
    ts = Vector{typeof(out[1][1])}(0)
    push!(ts,t)
    for (t,u,du) in out
      push!(ts,t)
      if typeof(u₀) <: AbstractArray
        push!(timeseries,u)
      else
        push!(timeseries,u[1])
      end
    end
    t = ts[end]
    u = timeseries[end]
  elseif alg ∈ SUNDIALS_ALGORITHMS
    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    # Needs robustness
    o[:Ts] = map(Float64,o[:Ts])
    o[:t] = map(Float64,o[:t])
    t = o[:t]; Ts = o[:Ts]
    initialize_backend(:Sundials)
    opts = buildOptions(o,SUNDIALS_OPTION_LIST,SUNDIALS_ALIASES,SUNDIALS_ALIASES_REVERSED)
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = vec(prob.f(reshape(u,sizeu),t)); 0)
    else
      f! = (t,u,du) -> (prob.f(reshape(du,sizeu),reshape(u,sizeu),t); u = vec(u); du=vec(du); 0)
    end
    ts = [t;Ts] # No adaptive
    vectimeseries = Sundials.cvode(f!,vec(u),ts)
    t = ts[end]
    if typeof(u₀)<:AbstractArray
      timeseries = GrowableArray(u₀;initvalue=false)
      for i=1:size(vectimeseries,1)
        push!(timeseries,reshape(vectimeseries[i,:]',sizeu))
      end
    else
      timeseries = vectimeseries
    end
    u = timeseries[end]
  end

  (atomloaded && progressbar) ? Main.Atom.progress(1) : nothing #Use Atom's progressbar if loaded

  if knownanalytic
    u_analytic = analytic(u₀,t)
    if save_timeseries
      timeseries_analytic = GrowableArray(analytic(u₀,ts[1]))
      for i in 2:size(timeseries,1)
        push!(timeseries_analytic,analytic(u₀,ts[i]))
      end
      timeseries = copy(timeseries)
      timeseries_analytic = copy(timeseries_analytic)
      return(ODESolution(u,u_analytic,timeseries=timeseries,t=ts,timeseries_analytic=timeseries_analytic))
    else
      return(ODESolution(u,u_analytic))
    end
  else #No known analytic
    if save_timeseries
      timeseries = copy(timeseries)
      return(ODESolution(u,timeseries=timeseries,t=ts))
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

function ode_determine_initΔt(u₀::AbstractArray,t,abstol,reltol,internalnorm,f,order)
  f₀ = similar(u₀); f₁ = similar(u₀); u₁ = similar(u₀)
  d₀ = norm(u₀./(abstol+u₀*reltol),internalnorm)
  f(f₀,u₀,t)
  d₁ = norm(f₀./(abstol+u₀*reltol),internalnorm)
  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    Δt₀ = 1//10^(6)
  else
    Δt₀ = (d₀/d₁)/100
  end
  @inbounds for i in eachindex(u₀)
     u₁[i] = u₀[i] + Δt₀*f₀[i]
  end
  f(f₁,u₁,t+Δt₀)
  d₂ = norm((f₁-f₀)./(abstol+u₀*reltol),internalnorm)/Δt₀
  if max(d₁,d₂)<=1//10^(15)
    Δt₁ = max(1//10^(6),Δt₀*1//10^(3))
  else
    Δt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order+1))
  end
  Δt = min(100Δt₀,Δt₁)
end

function ode_determine_initΔt(u₀::Number,t,abstol,reltol,internalnorm,f,order)
  d₀ = norm(u₀./(abstol+u₀*reltol),internalnorm)
  f₀ =f(u₀,t)
  d₁ = norm(f₀./(abstol+u₀*reltol),internalnorm)
  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    Δt₀ = 1//10^(6)
  else
    Δt₀ = (d₀/d₁)/100
  end
  u₁ = u₀ + Δt₀*f₀
  f₁ = f(u₁,t+Δt₀)
  d₂ = norm((f₁-f₀)./(abstol+u₀*reltol),internalnorm)/Δt₀
  if max(d₁,d₂)<=1//10^(15)
    Δt₁ = max(1//10^(6),Δt₀*1//10^(3))
  else
    Δt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order+1))
  end
  Δt = min(100Δt₀,Δt₁)
end
