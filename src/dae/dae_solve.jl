"""
`solve(prob::DAEProblem,tspan)`

Solves the DAE as defined by prob on the time interval tspan. If not given, tspan defaults to [0,1].

### Keyword Arguments

* `Δt`: Sets the initial stepsize. Defaults to an automatic choice.
* `save_timeseries`: Saves the result at every timeseries_steps steps. Default is true.
* `timeseries_steps`: Denotes how many steps between saving a value for the timeseries. Defaults to 1.
* `adaptive` - Turns on adaptive timestepping for appropriate methods. Default is false.
* `γ` - The risk-factor γ in the q equation for adaptive timestepping. Default is 2.
* `qmax` - Defines the maximum value possible for the adaptive q. Default is 1.125.
* `ablstol` - Absolute tolerance in adaptive timestepping. Defaults to 1e-3.
* `reltol` - Relative tolerance in adaptive timestepping. Defaults to 1e-6.
* `maxiters` - Maximum number of iterations before stopping. Defaults to 1e9.
* `Δtmax` - Maximum Δt for adaptive timestepping. Defaults to half the timespan.
* `Δtmin` - Minimum Δt for adaptive timestepping. Defaults to 1e-10.
* `internalnorm` - The norm for which error estimates are calculated. Default is 2.
* `progressbar` - Turns on/off the Juno progressbar. Defualt is false.
* `progress_steps` - Numbers of steps between updates of the progress bar. Default is 1000.
* `alg`: String which defines the solver algorithm. Defult is "idasol". Possibilities are:
  - `idasol`: The DAE solver from Sundials

"""
function solve(prob::AbstractDAEProblem,tspan::AbstractArray=[0,1];Δt::Number=0.0,save_timeseries::Bool = true,
              timeseries_steps::Int = 1,alg=nothing,adaptive=false,γ=2.0,alg_hint=nothing,
              abstol=1e-3,reltol=1e-6,qmax=1.125,maxiters::Int = round(Int,1e9),
              Δtmax=nothing,Δtmin=nothing,progress_steps=1000,internalnorm=2, saveat=[],
              progressbar=false,tType=typeof(Δt))

  if tspan[end]-tspan[1]<0
    tspan = vec(tspan)
    error("final time must be greater than starting time. Aborting.")
  end
  atomloaded = isdefined(Main,:Atom)
  t = tspan[1]
  Ts = tspan[2:end]
  @unpack u₀,du₀,knownanalytic,analytic,numvars,isinplace = prob
  uType = typeof(u₀)
  uEltype = eltype(u₀)
  rateType = typeof(du₀)

  u = copy(u₀)
  du= copy(du₀)
  ks = Vector{uType}(0)

  if alg == nothing
    alg = plan_dae(alg_hint,abstol,reltol)
  end


  if alg == :idasol
    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    u = map(Float64,u) # Needs Float64
    # Needs robustness
    Ts = map(Float64,Ts)
    t = map(Float64,t)
    saveat = [float(x) for x in saveat]
    initialize_backend(:Sundials)
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du,out) -> (du[:] = vec(prob.f(t,reshape(u,sizeu),reshape(du,sizeu))); 0)
    else
      f! = (t,u,du,out) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu),reshape(out,sizeu)); u = vec(u); du=vec(du); out=vec(out); 0)
    end
    ts = [t;Ts]

    vectimeseries,vectimeseries_du = Sundials.idasol(f!,u,du,ts)
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
    t = ts[end]
    u = timeseries[end]
  end

  if knownanalytic
    u_analytic,du_analytic = analytic(t,u₀,du₀)
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      u_tmp,du_tmp = analytic(ts[i],u₀,du₀)
      push!(timeseries_analytic,u_tmp)
      push!(timeseries_du_analytic,du_tmp)
    end
    return(DAESolution(u,u_analytic,du_analytic,prob,alg,timeseries=timeseries,t=ts,
            timeseries_analytic=timeseries_analytic,
            timeseries_du_analytic=timeseries_du_analytic,
            k=ks,saveat=saveat))
  else
    return(DAESolution(u,du,prob,alg,timeseries=timeseries,t=ts,k=ks,saveat=saveat))
  end
end


function plan_dae(alg_hint,abstol,reltol)
  :idasol
end
