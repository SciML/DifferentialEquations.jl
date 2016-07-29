"""
solve(prob::SDEProblem,tspan=[0,1];Δt=0)

Solves the SDE as defined by prob with initial Δt on the time interval tspan.
If not given, tspan defaults to [0,1]. If

### Keyword Arguments

* `Δt`: Sets the initial stepsize. Defaults to an automatic choice.
* `save_timeseries`: Saves the result at every timeseries_steps steps. Default is true.
* `timeseries_steps`: Denotes how many steps between saving a value for the timeseries. Defaults to 1.
* `adaptive` - Turns on adaptive timestepping for appropriate methods. Default is false.
* `γ` - The risk-factor γ in the q equation for adaptive timestepping. Default is 2.
* `qmax` - Defines the maximum value possible for the adaptive q. Default is 1.125.
* `δ` - The weight-factor in the error estimate. Default is 1/6.
* `ablstol` - Absolute tolerance in adaptive timestepping. Defaults to 1e-3.
* `reltol` - Relative tolerance in adaptive timestepping. Defaults to 1e-6.
* `maxiters` - Maximum number of iterations before stopping. Defaults to 1e9.
* `Δtmax` - Maximum Δt for adaptive timestepping. Defaults to half the timespan.
* `Δtmin` - Minimum Δt for adaptive timestepping. Defaults to 1e-10.
* `internalnorm` - The norm for which error estimates are calculated. Default is 2.
* `progressbar` - Turns on/off the Juno progressbar. Defualt is false.
* `progress_steps` - Numbers of steps between updates of the progress bar. Default is 1000.
* `discard_length` - Size at which to discard future information in adaptive. Default is 1e-15.
* `tableau`: The tableau for an `:SRA` or `:SRI` algorithm. Defaults to SRIW1 or SRA1.
* `adaptivealg`: The adaptive timestepping algorithm. Default is `:RSwm3`.
* `alg`: String which defines the solver algorithm. Defult is "SRIW1Optimized". Possibilities are:
    - `:EM`- The Euler-Maruyama method.
    - `:RKMil` - An explicit Runge-Kutta discretization of the strong Order 1.0 Milstein method.
    - `:SRA` - The strong Order 2.0 methods for additive SDEs due to Rossler. Not yet implemented.
        Default tableau is for SRA1.
    - `:SRI` - The strong Order 1.5 methods for diagonal/scalar SDEs due to Rossler.
  Default tableau is for SRIW1.
    - `:SRIW1Optimized` - An optimized version of SRIW1. Strong Order 1.5.
    - `:SRA1Optimized` - An optimized version of SRIA1. Strong Order 2.0.
    - `:SRAVectorized` - A vectorized implementation of SRA algorithms. Requires 1-dimensional problem.
    - `:SRIVectorized` - A vectorized implementation of SRI algorithms. Requires 1-dimensional problem.
"""
function solve(prob::SDEProblem,tspan::AbstractArray=[0,1];Δt::Number=0,save_timeseries::Bool = true,
              timeseries_steps::Int = 1,alg::Symbol=:SRIW1Optimzed,adaptive=false,γ=2.0,
              abstol=1e-3,reltol=1e-6,qmax=1.125,δ=1/6,maxiters::Int = round(Int,1e9),
              Δtmax=nothing,Δtmin=nothing,progress_steps=1000,internalnorm=2,
              discard_length=1e-15,adaptivealg::Symbol=:RSwM3,progressbar=false,tType=typeof(Δt),tableau = nothing)

  @unpack prob: f,σ,u₀,knownsol,sol, numvars, sizeu

  tspan = vec(tspan)
  if tspan[2]-tspan[1]<0 || length(tspan)>2
    error("tspan must be two numbers and final time must be greater than starting time. Aborting.")
  end

  if adaptive && alg ∈ SDE_ADAPTIVEALGORITHMS
    tType = Float64
    initialize_backend(:DataStructures)
    if adaptivealg == :RSwM3
      initialize_backend(:ResettableStacks)
    end
  end

  if Δt == 0.0
    if alg==:Euler
      order = 0.5
    elseif alg==:RKMil
      order = 1.0
    else
      order = 1.5
    end
    Δt = sde_determine_initΔt(u₀,float(tspan[1]),abstol,reltol,internalnorm,f,σ,order)
    tType=typeof(Δt)
  end

  if Δtmax == nothing
    Δtmax = (tspan[2]-tspan[1])/2
  end
  if Δtmin == nothing
    Δtmin = tType(1e-10)
  end


  uType = typeof(u₀)

  T = tType(tspan[2])
  t = tType(tspan[1])
  u = u₀
  if numvars == 1
    W = 0.0
    Z = 0.0
  else
    W = zeros(sizeu)
    Z = zeros(sizeu)
  end



  timeseries = GrowableArray(u)
  ts = Vector{tType}(0)
  Ws = GrowableArray(W)
  push!(ts,t)

  #PreProcess
  if alg==:SRA && tableau == nothing
    tableau = constructSRA1()
  elseif alg==:SRI && tableau == nothing
    tableau = constructSRIW1()
  end

  if numvars == 1
    rands = ChunkedArray(randn)
  else
    rands = ChunkedArray(randn,u)
  end

  sqΔt = sqrt(Δt)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  maxStackSize = 0
  #EEst = 0

  # This part is a mess. Needs cleaning.
  if alg==:EM
    u,t,W,timeseries,ts,Ws = sde_eulermaruyama(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,progressbar,rands,sqΔt,W)
  elseif alg==:RKMil
    u,t,W,timeseries,ts,Ws = sde_rkmil(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,progressbar,rands,sqΔt,W)
  elseif alg==:SRI
    u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2 =            sde_sri(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,rands,sqΔt,W,Z,tableau)
  elseif alg==:SRIW1Optimized
    u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2 = sde_sriw1optimized(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,rands,sqΔt,W,Z,tableau)
  elseif alg==:SRIVectorized
    u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2 =  sde_srivectorized(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,rands,sqΔt,W,Z,tableau)
  elseif alg==:SRAVectorized
    u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2 =  sde_sravectorized(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,rands,sqΔt,W,Z,tableau)
  elseif alg==:SRA1Optimized || alg==:SRA # Need a devectorized form
    u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2 =  sde_sra1optimized(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,rands,sqΔt,W,Z,tableau)
  end

  (atomloaded && progressbar) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded

  if knownsol
    uTrue = sol(u₀,t,W)
    if save_timeseries
      sols = GrowableArray(sol(u₀,ts[1],Ws[1]))
      for i in 2:size(Ws,1)
        push!(sols,sol(u₀,ts[i],Ws[i]))
      end
      Ws = copy(Ws)
      timeseries = copy(timeseries)
      sols = copy(sols)
      return(SDESolution(u,uTrue,W=W,timeseries=timeseries,ts=ts,Ws=Ws,sols=sols,maxStackSize=maxStackSize))
    else
      return(SDESolution(u,uTrue,W=W,maxStackSize=maxStackSize))
    end
  else #No known sol
    if save_timeseries
      timeseries = copy(timeseries)
      return(SDESolution(u,timeseries=timeseries,W=W,ts=ts,maxStackSize=maxStackSize))
    else
      return(SDESolution(u,W=W,maxStackSize=maxStackSize))
    end
  end
end

const SDE_ADAPTIVEALGORITHMS = Set([:SRI,:SRIW1Optimized,:SRIVectorized,:SRAVectorized,:SRA1Optimized,:SRA])

function sde_determine_initΔt(u₀,t,abstol,reltol,internalnorm,f,σ,order)
  d₀ = norm(u₀./(abstol+u₀*reltol),2)
  f₀ = f(u₀,t)
  σ₀ = 3σ(u₀,t)

  d₁ = norm(max(abs(f₀+σ₀),abs(f₀-σ₀))./(abstol+u₀*reltol),2)
  if d₀ < 1e-5 || d₁ < 1e-5
    Δt₀ = 1e-6
  else
    Δt₀ = 0.01*(d₀/d₁)
  end
  u₁ = u₀ + Δt₀*f₀
  f₁ = f(u₁,t+Δt₀)
  σ₁ = 3σ(u₁,t+Δt₀)
  ΔσMax = max(abs(σ₀-σ₁),abs(σ₀+σ₁))
  d₂ = norm(max(abs(f₁-f₀+ΔσMax),abs(f₁-f₀-ΔσMax))./(abstol+u₀*reltol),2)/Δt₀
  if max(d₁,d₂)<=1e-15
    Δt₁ = max(1e-6,Δt₀*1e-3)
  else
    if !isdefined(Main,:order)
      order = 1 #Convervative choice
    end
    Δt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order+.5))
  end
  Δt = min(100Δt₀,Δt₁)
end
