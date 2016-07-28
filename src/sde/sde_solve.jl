"""
solve(prob::SDEProblem,tspan=[0,1];Δt=0)

Solves the SDE as defined by prob with initial Δt on the time interval tspan.
If not given, tspan defaults to [0,1]. If

### Keyword Arguments

* fullSave: Saves the result at every saveSteps steps. Default is false.
saveSteps: If fullSave is true, then the output is saved every saveSteps steps.
* alg: String which defines the solver algorithm. Defult is "SRI". Possibilities are:
  * "EM"- The Euler-Maruyama method.
  * "RKMil" - An explicit Runge-Kutta discretization of the strong Order 1.0 Milstein method.
  * "SRA" - The strong Order 1.5 method for additive SDEs due to Rossler.
  * "SRI" - The strong Order 1.5 method for diagonal/scalar SDEs due to Rossler. Most efficient.
"""
function solve(prob::SDEProblem,tspan::AbstractArray=[0,1];Δt::Number=0,fullSave::Bool = false,
              saveSteps::Int = 1,alg::Symbol=:SRI,adaptive=false,γ=2.0,
              abstol=1e-3,reltol=1e-2,qmax=1.125,δ=1/6,maxiters::Int = round(Int,1e15),
              Δtmax=nothing,Δtmin=nothing,progressSteps=1000,internalNorm=2,
              discardLength=1e-15,adaptivealg::Symbol=:RSwM3,progressBar=false,tType=typeof(Δt),tableau = constructSRIW1())

  @unpack prob: f,σ,u₀,knownSol,sol, numVars, sizeu

  tspan = vec(tspan)
  if tspan[2]-tspan[1]<0 || length(tspan)>2
    error("tspan must be two numbers and final time must be greater than starting time. Aborting.")
  end

  if adaptive && alg ∈ SDE_ADAPTIVEALGORITHMS
    tType = Float64
  end

  if Δt == 0.0
    if alg==:Euler
      order = 0.5
    elseif alg==:RKMil
      order = 1.0
    else
      order = 1.5
    end
    Δt = sde_determine_initΔt(u₀,float(tspan[1]),abstol,reltol,internalNorm,f,σ,order)
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
  if numVars == 1
    W = 0.0
    Z = 0.0
  else
    W = zeros(sizeu)
    Z = zeros(sizeu)
  end



  uFull = GrowableArray(u)
  tFull = Vector{tType}(0)
  WFull = GrowableArray(W)
  push!(tFull,t)

  #PreProcess

  if alg==:SRA
    SRA = constructSRA1()
  end

  if numVars == 1
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

  if alg==:EM
    u,t,W,uFull,tFull,WFull = sde_eulermaruyama(f,σ,u,t,Δt,T,iter,maxiters,uFull,WFull,tFull,saveSteps,fullSave,adaptive,progressBar,rands,sqΔt,W)
  elseif alg==:RKMil
    u,t,W,uFull,tFull,WFull = sde_rkmil(f,σ,u,t,Δt,T,iter,maxiters,uFull,WFull,tFull,saveSteps,fullSave,adaptive,progressBar,rands,sqΔt,W)
  elseif alg==:SRI
    u,t,W,uFull,tFull,WFull,maxStackSize,maxStackSize2 = sde_sri(f,σ,u,t,Δt,T,iter,maxiters,uFull,WFull,tFull,saveSteps,fullSave,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,numVars,sizeu,discardLength,progressBar,rands,sqΔt,W,Z,tableau,)
  elseif alg==:SRIW1Optimized
    u,t,W,uFull,tFull,WFull,maxStackSize,maxStackSize2 = sde_sriw1optimized(f,σ,u,t,Δt,T,iter,maxiters,uFull,WFull,tFull,saveSteps,fullSave,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,numVars,sizeu,discardLength,progressBar,rands,sqΔt,W,Z,tableau)
  elseif alg==:SRIVectorized
    u,t,W,uFull,tFull,WFull,maxStackSize,maxStackSize2 = sde_srivectorized(f,σ,u,t,Δt,T,iter,maxiters,uFull,WFull,tFull,saveSteps,fullSave,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,numVars,sizeu,discardLength,progressBar,rands,sqΔt,W,Z,tableau)
  elseif alg==:SRAVectorized
    u,t,W,uFull,tFull,WFull,maxStackSize,maxStackSize2 = sde_sravectorized(f,σ,u,t,Δt,T,iter,maxiters,uFull,WFull,tFull,saveSteps,fullSave,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,numVars,sizeu,discardLength,progressBar,rands,sqΔt,W,Z,tableau)
  elseif alg==:SRA1Optimized || alg==:SRA # Need a devectorized form
    u,t,W,uFull,tFull,WFull,maxStackSize,maxStackSize2 = sde_sra1optimized(f,σ,u,t,Δt,T,iter,maxiters,uFull,WFull,tFull,saveSteps,fullSave,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,numVars,sizeu,discardLength,progressBar,rands,sqΔt,W,Z,tableau)
  end
  (atomLoaded && progressBar) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded

  if knownSol
    uTrue = sol(u₀,t,W)
    if fullSave
      solFull = GrowableArray(sol(u₀,tFull[1],WFull[1]))
      for i in 2:size(WFull,1)
        push!(solFull,sol(u₀,tFull[i],WFull[i]))
      end
      WFull = copy(WFull)
      uFull = copy(uFull)
      solFull = copy(solFull)
      return(SDESolution(u,uTrue,W=W,uFull=uFull,tFull=tFull,WFull=WFull,solFull=solFull,maxStackSize=maxStackSize))
    else
      return(SDESolution(u,uTrue,W=W,maxStackSize=maxStackSize))
    end
  else #No known sol
    if fullSave
      uFull = copy(uFull)
      return(SDESolution(u,uFull=uFull,W=W,tFull=tFull,maxStackSize=maxStackSize))
    else
      return(SDESolution(u,W=W,maxStackSize=maxStackSize))
    end
  end
end

const SDE_ADAPTIVEALGORITHMS = Set([:SRI,:SRIW1Optimized,:SRIVectorized,:SRAVectorized,:SRA1Optimized,:SRA])

function sde_determine_initΔt(u₀,t,abstol,reltol,internalNorm,f,σ,order)
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
