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
* `abstol` - Absolute tolerance in adaptive timestepping. Defaults to 1e-3.
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
  if tspan[end]-tspan[1]<0
    tspan = vec(tspan)
    error("final time must be greater than starting time. Aborting.")
  end
  atomloaded = isdefined(Main,:Atom)
  o = KW(kwargs)
  o[:t] = tspan[1]
  o[:Ts] = tspan[2:end]
  @unpack u₀,knownanalytic,analytic,numvars,isinplace = prob

  command_opts = copy(DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS)
  for (k,v) in o
    command_opts[k]=v
  end
  # Get the control variables
  @materialize save_timeseries, progressbar = command_opts

  if command_opts[:callback] == nothing
    callback = ODE_DEFAULT_CALLBACK
    custom_callback = false
  else
    callback = command_opts[:callback]
    custom_callback = true
  end

  u = copy(u₀)
  ks = Vector{uType}(0)
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
      @unpack order,adaptiveorder = o[:tableau]
    end
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = prob.f(t,u))
    else
      f! = prob.f
    end
    if Δt==0
      Δt = ode_determine_initΔt(u₀,float(tspan[1]),o[:abstol],o[:reltol],o[:internalnorm],f!,order)
    end
    if alg ∉ DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS
      o[:adaptive] = false
    else
      if o[:adaptive] == true
        Δt = 1.0*Δt # Convert to float in a way that keeps units
      end
    end
    if alg ∈ DIFFERENTIALEQUATIONSJL_IMPLICITALGS
      initialize_backend(:NLsolve)
    end

    if o[:tType] == nothing # if tType is not specified, grab it from Δt which defaults to 0.0 => Float64
      tType=typeof(Δt)
    else
      tType = o[:tType]
    end

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
    saveat = tType[convert(tType,x) for x in o[:saveat]]


    ### Algorithm-specific defaults ###

    if o[:qmin] == nothing # Use default qmin
      if alg == :DP5 || alg == :DP5Vectorized || alg == :DP5Threaded
        qmin = 0.2
      elseif alg == :DP8 || alg == :DP8Vectorized
        qmin = 0.333
      else
        qmin = 0.2
      end
    else
      qmin = o[:qmin]
    end
    if o[:qmax] == nothing # Use default qmax
      if alg == :DP5 || alg == :DP5Vectorized || alg == :DP5Threaded
        qmax = 10.0
      elseif alg == :DP8 || alg == :DP8Vectorized
        qmax = 6.0
      else
        qmax = 10.0
      end
    else
      qmax = o[:qmax]
    end
    if o[:β] == nothing # Use default β
      if alg == :DP5 || alg == :DP5Vectorized || alg == :DP5Threaded
        β = 0.04 # More than Hairer's suggestion
      elseif alg == :DP8 || alg == :DP8Vectorized
        β = 0.07 # More than Hairer's suggestion
      else
        β = 0.4 / order
      end
    else
      β = o[:β]
    end
    if o[:expo1] == nothing # Use default expo1
      if alg == :DP5 || alg == :DP5Vectorized || alg == :DP5Threaded
        expo1 = 1/order - .75β
      elseif alg == :DP8 || alg == :DP8Vectorized
        expo1 = 1/order - .2β
      else
        expo1 = .7/order
      end
    else
      expo1 = o[:expo1]
    end
    fsal = false
    if alg ∈ DIFFERENTIALEQUATIONSJL_FSALALGS
      fsal = true
    elseif alg == :ExplicitRK
      @unpack fsal = o[:tableau]
    end

    Ts = map(tType,o[:Ts])
    t = tType(o[:t])

    o[:abstol] = convert(uEltype,o[:abstol])

    if uType <: AbstractArray
      uEltypeNoUnits = eltype(u./u)
    else
      uEltypeNoUnits = typeof(u./u)
    end
    rateType = typeof(u/zero(t)) ## Can be different if united

    if alg ∈ DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS
      ksEltype = Vector{rateType} # Store more ks for the special algs
    else
      ksEltype = rateType # Makes simple_dense
    end

    @materialize maxiters,timeseries_steps,save_timeseries,adaptive,progress_steps,abstol,reltol,γ,Δtmax,Δtmin,internalnorm,tableau,autodiff,timechoicealg,qoldinit,dense = o
    #println(@code_typed ode_solve(ODEIntegrator{alg,uType,uEltype,ndims(u)+1,tType,uEltypeNoUnits,rateType,ksEltype}(f!,u,t,Δt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,Δtmax,Δtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β,expo1,timechoicealg,qoldinit,normfactor,fsal,dense,saveat,alg,callback,custom_callback)))

    u,t,timeseries,ts,ks = ode_solve(ODEIntegrator{alg,uType,uEltype,ndims(u)+1,tType,uEltypeNoUnits,rateType,ksEltype}(f!,u,t,Δt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,Δtmax,Δtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β,expo1,timechoicealg,qoldinit,normfactor,fsal,dense,saveat,alg,callback,custom_callback))

    if ts[end] != t
      push!(ts,t)
      push!(timeseries,u)
    end
  elseif alg ∈ ODEINTERFACE_ALGORITHMS
    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    o[:Ts] = float(o[:Ts])
    o[:t] = float(o[:t])
    t = o[:t]; Ts = o[:Ts]
    saveat = [float(x) for x in command_opts[:saveat]]
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = vec(prob.f(t,reshape(u,sizeu))); nothing)
    else
      f! = (t,u,du) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu)); u = vec(u); du=vec(du); nothing)
    end
    initialize_backend(:ODEInterface)
    o[:RHS_CALLMODE] = ODEInterface.RHS_CALL_INSITU
    dict = buildOptions(o,ODEINTERFACE_OPTION_LIST,ODEINTERFACE_ALIASES,ODEINTERFACE_ALIASES_REVERSED)
    opts = ODEInterface.OptionsODE([Pair(ODEINTERFACE_STRINGS[k],v) for (k,v) in dict]...) #Convert to the strings
    du = similar(u)
    if alg==:dopri5
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dopri5,f!,[t;Ts],vec(u),opts)
    elseif alg==:dop853
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dop853,f!,[t;Ts],vec(u),opts)
    elseif alg==:odex
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.odex,f!,[t;Ts],vec(u),opts)
    elseif alg==:seulex
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.seulex,f!,[t;Ts],vec(u),opts)
    elseif alg==:radau
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau,f!,[t;Ts],vec(u),opts)
    elseif alg==:radau5
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau5,f!,[t;Ts],vec(u),opts)
    end
    t = ts[end]
    if typeof(u₀)<:AbstractArray
      timeseries = Vector{uType}(0)
      for i=1:size(vectimeseries,1)
        push!(timeseries,reshape(view(vectimeseries,i,:,)',sizeu))
      end
    else
      timeseries = vec(vectimeseries)
    end
    u = timeseries[end]

  elseif alg ∈ ODEJL_ALGORITHMS
    if typeof(u) <: Number
      u = [u]
    end
    # Needs robustness
    o[:Ts] = float(o[:Ts])
    o[:t] = float(o[:t])
    t = o[:t]; Ts = o[:Ts]
    o[:T] = Ts[end]
    saveat = [float(x) for x in command_opts[:saveat]]
    initialize_backend(:ODEJL)
    opts = buildOptions(o,ODEJL_OPTION_LIST,ODEJL_ALIASES,ODEJL_ALIASES_REVERSED)
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = prob.f(t,u))
    else
      f! = prob.f
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
    out = ODE.solve(ode;solver=solver,opts...)
    timeseries = out.y
    ts = out.t
    if length(out.y[1])==1
      tmp = Vector{eltype(out.y[1])}(length(out.y))
      for i in 1:length(out.y)
        tmp[i] = out.y[i][1]
      end
      timeseries = tmp
    end
    t = ts[end]
    u = timeseries[end]
  elseif alg ∈ SUNDIALS_ALGORITHMS
    if alg == :cvode_BDF
      integrator = :BDF
    elseif alg ==  :cvode_Adams
      integrator = :Adams
    end

    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    u = map(Float64,u) # Needs Float64
    # Needs robustness
    o[:Ts] = map(Float64,o[:Ts])
    o[:t] = map(Float64,o[:t])
    t = o[:t]; Ts = o[:Ts];
    saveat = [float(x) for x in command_opts[:saveat]]
    initialize_backend(:Sundials)
    opts = buildOptions(o,SUNDIALS_OPTION_LIST,SUNDIALS_ALIASES,SUNDIALS_ALIASES_REVERSED)
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = vec(prob.f(t,reshape(u,sizeu))); 0)
    else
      f! = (t,u,du) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu)); u = vec(u); du=vec(du); 0)
    end
    ts = [t;Ts]
    @materialize abstol, reltol = command_opts
    if command_opts[:adaptive]
      ts, vectimeseries = Sundials.cvode_fulloutput(f!,vec(u),ts;integrator=integrator,abstol=float(abstol),reltol=float(reltol))
      timeseries = Vector{uType}(0)
      if typeof(u₀)<:AbstractArray
        for i=1:size(vectimeseries,1)
          push!(timeseries,reshape(vectimeseries[i],sizeu))
        end
      else
        for i=1:size(vectimeseries,1)
          push!(timeseries,vectimeseries[i][1])
        end
      end
    else
      Δt = command_opts[:Δt]
      ts = float(collect(t:Δt:Ts[end]))
      if length(Ts)>1
        ts = float([ts;Ts[1:end-1]])
        sort(ts)
      end
      vectimeseries = Sundials.cvode(f!,vec(u),ts,integrator=integrator)
      timeseries = Vector{uType}(0)
      if typeof(u₀)<:AbstractArray
        for i=1:size(vectimeseries,1)
          push!(timeseries,reshape(view(vectimeseries,i,:),sizeu))
        end
      else
        for i=1:size(vectimeseries,1)
          push!(timeseries,vectimeseries[i])
        end
      end
    end
    t = ts[end]
    u = timeseries[end]
  end

  (atomloaded && progressbar) ? Main.Atom.progress(1) : nothing #Use Atom's progressbar if loaded

  if knownanalytic
    u_analytic = analytic(t,u₀)
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      push!(timeseries_analytic,analytic(ts[i],u₀))
    end
    return(ODESolution(u,u_analytic,prob,alg,timeseries=timeseries,t=ts,timeseries_analytic=timeseries_analytic,k=ks,saveat=saveat))
  else
    return(ODESolution(u,prob,alg,timeseries=timeseries,t=ts,k=ks,saveat=saveat))
  end
  #=

  if knownanalytic
    u_analytic = analytic(t,u₀)
    if save_timeseries
      timeseries_analytic = Vector{uType}(0)
      for i in 1:size(timeseries,1)
        push!(timeseries_analytic,analytic(ts[i],u₀))
      end
      return(ODESolution(u,u_analytic,prob,alg,timeseries=timeseries,t=ts,timeseries_analytic=timeseries_analytic,k=ks,saveat=saveat))
    else
      return(ODESolution(u,u_analytic,prob,alg))
    end
  else #No known analytic
    if save_timeseries
      return(ODESolution(u,prob,alg,timeseries=timeseries,t=ts,k=ks,saveat=saveat))
    else
      return(ODESolution(u,prob,alg))
    end
  end
  =#
end

function buildOptions(o,optionlist,aliases,aliases_reversed)
  dict1 = Dict{Symbol,Any}([Pair(k,o[k]) for k in (keys(o) ∩ optionlist)])
  dict2 = Dict([Pair(aliases_reversed[k],o[k]) for k in (keys(o) ∩ values(aliases))])
  merge(dict1,dict2)
end

function ode_determine_initΔt(u₀::AbstractArray,t,abstol,reltol,internalnorm,f,order)
  f₀ = similar(u₀); f₁ = similar(u₀); u₁ = similar(u₀)
  d₀ = norm(u₀./(abstol+u₀*reltol),internalnorm)
  f(t,u₀,f₀)
  d₁ = norm(f₀./(abstol+u₀*reltol),internalnorm)
  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    Δt₀ = 1//10^(6)
  else
    Δt₀ = (d₀/d₁)/100
  end
  @inbounds for i in eachindex(u₀)
     u₁[i] = u₀[i] + Δt₀*f₀[i]
  end
  f(t+Δt₀,u₁,f₁)
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
  f₀ =f(t,u₀)
  d₁ = norm(f₀./(abstol+u₀*reltol),internalnorm)
  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    Δt₀ = 1//10^(6)
  else
    Δt₀ = (d₀/d₁)/100
  end
  u₁ = u₀ + Δt₀*f₀
  f₁ = f(t+Δt₀,u₁)
  d₂ = norm((f₁-f₀)./(abstol+u₀*reltol),internalnorm)/Δt₀
  if max(d₁,d₂)<=1//10^(15)
    Δt₁ = max(1//10^(6),Δt₀*1//10^(3))
  else
    Δt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order+1))
  end
  Δt = min(100Δt₀,Δt₁)
end
