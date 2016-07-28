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
    if Δt==0
      Δt = determine_initΔt(u₀,abstol,reltol,internalNorm,f,order)
    end
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
    elseif alg==:dopri853
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

const DIFFERENTIALEQUATIONSJL_ALGORITHMS = Set([:Euler,:Midpoint,:RK4,:ExplicitRK,:ImplicitEuler,:Trapezoid,:Rosenbrock32])
const ODEINTERFACE_ALGORITHMS = Set([:dopri5,:dopri853,:odex,:radau5,:radau,:seulex])
const ODEJL_ALGORITHMS = Set([:ode23,:ode45,:ode78,:ode23s,:ode1,:ode2_midpoint,:ode2_heun,:ode4,:ode45_fe])

const DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS = Dict(:Δt => 0,
                                 :fullSave => false,
                                 :saveSteps => 1,
                                 :tableau => DifferentialEquations.DEFAULT_TABLEAU,
                                 :adaptive => true,
                                 :γ=>2.0,
                                 :abstol=>1//10^8,
                                 :reltol=>1//10^6,
                                 :qmax=>4,
                                 :maxiters => round(Int,1e9),
                                 :Δtmax=>nothing,
                                 :Δtmin=>nothing,
                                 :internalNorm => 2,
                                 :progressBar=>false,
                                 :progressSteps=>1000)

const ODEJL_OPTION_LIST = Set([:tout,:tstop,:reltol,:abstol,:minstep,:maxstep,:initstep,:norm,:maxiters,:isoutofdomain])
const ODEJL_ALIASES = Dict{Symbol,Symbol}(:minstep=>:Δtmin,:maxstep=>:Δtmax,:initstep=>:Δt,:tstop=>:T)
const ODEJL_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEJL_ALIASES])
const ODEINTERFACE_OPTION_LIST = Set([:RTOL,:ATOL,:OUTPUTFCN,:OUTPUTMODE,
                                :MAXSTEPS,:STEST,:EPS,:RHO,:SSMINSEL,
                                :SSMAXSEL,:SSBETA,:MAXSS,:INITIALSS,
                                :MAXEXCOLUMN,:STEPSIZESEQUENCE,:MAXSTABCHECKS,
                                :MAXSTABCHECKLINE,:DENSEOUTPUTWOEE,:INTERPOLDEGRE,
                                :SSREDUCTION,:SSSELECTPAR1,:SSSELECTPAR2,
                                :ORDERDECFRAC,:ORDERINCFRAC,:OPT_RHO,:OPT_RHO2,
                                :RHSAUTONOMOUS,:M1,:M2,:LAMBDADENSE,:TRANSJTOH,
                                :STEPSIZESEQUENCE,:JACRECOMPFACTOR,:MASSMATRIX,
                                :JACOBIMATRIX,:JACOBIBANDSSTRUCT,:WORKFORRHS,
                                :WORKFORJAC,:WORKFORDEC,:WORKFORSOL,:MAXNEWTONITER,
                                :NEWTONSTARTZERO,:NEWTONSTOPCRIT,:DIMFIND1VAR,
                                :MAXSTAGES,:MINSTAGES,:INITSTAGES,:STEPSIZESTRATEGY,
                                :FREEZESSLEFT,:FREEZESSRIGHT,:ORDERDECFACTOR,
                                :ORDERINCFACTOR,:ORDERDECCSTEPFAC1,:ORDERDECSTEPFAC2
                                ])

const ODEINTERFACE_ALIASES = Dict{Symbol,Symbol}(:RTOL=>:reltol,
                                                 :ATOL=>:abstol,
                                                 :MAXSTEPS=> :maxiters,
                                                 :MAXSS=>:Δtmax,
                                                 :INITIALSS=>:Δt,
                                                 #:SSMINSEL=>:qmin,
                                                 :SSMAXSEL=>:qmax)
const ODEINTERFACE_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEINTERFACE_ALIASES])

const DIFFERENTIALEQUATIONSJL_ORDERS = Dict{Symbol,Int}(:Euler=>1,
                                                        :Midpoint=>2,
                                                        :RK4=>4,
                                                        :ExplicitRK=>4, #Gets overwritten
                                                        :ImplicitEuler=>1,
                                                        :Trapezoid=>2,
                                                        :Rosenbrock32=>3)
const DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS = Set([:ExplicitRK,:Rosenbrock32])
const ODEINTERFACE_STRINGS = Dict{Symbol,String}(
  :LOGIO            => "logio",
  :LOGLEVEL         => "loglevel",
  :RHS_CALLMODE     => "RightHandSideCallMode",

  :RTOL             => "RelTol",
  :ATOL             => "AbsTol",
  :MAXSTEPS         => "MaxNumberOfSteps",
  :EPS              => "eps",

  :OUTPUTFCN        => "OutputFcn",
  :OUTPUTMODE       => "OutputFcnMode",

  :STEST            => "StiffTestAfterStep",
  :RHO              => "rho",
  :SSMINSEL         => "StepSizeMinSelection",
  :SSMAXSEL         => "StepSizeMaxSelection",
  :SSBETA           => "StepSizeBeta",
  :MAXSS            => "MaxStep",
  :INITIALSS        => "InitialStep",


  :MAXEXCOLUMN      => "MaxExtrapolationColumn",
  :MAXSTABCHECKS    => "MaxNumberOfStabilityChecks",
  :MAXSTABCHECKLINE => "MaxLineForStabilityCheck",
  :INTERPOLDEGREE   => "DegreeOfInterpolation",
  :ORDERDECFRAC     => "OrderDecreaseFraction",
  :ORDERINCFRAC     => "OrderIncreaseFraction",
  :STEPSIZESEQUENCE => "StepSizeSequence",
  :SSREDUCTION      => "StepSizeReduction",
  :SSSELECTPAR1     => "StepSizeSelectionParam1",
  :SSSELECTPAR2     => "StepSizeSelectionParam2",
  :RHO2             => "rho2",
  :DENSEOUTPUTWOEE  => "DeactivateErrorEstInDenseOutput",

  :TRANSJTOH        => "TransfromJACtoHess",
  :MAXNEWTONITER    => "MaxNewtonIterations",
  :NEWTONSTARTZERO  => "StartNewtonWithZeros",
  :DIMOFIND1VAR     => "DimensionOfIndex1Vars",
  :DIMOFIND2VAR     => "DimensionOfIndex2Vars",
  :DIMOFIND3VAR     => "DimensionOfIndex3Vars",
  :STEPSIZESTRATEGY => "StepSizeStrategy",
  :M1               => "M1",
  :M2               => "M2",
  :JACRECOMPFACTOR  => "RecomputeJACFactor",
  :NEWTONSTOPCRIT   => "NewtonStopCriterion",
  :FREEZESSLEFT     => "FreezeStepSizeLeftBound",
  :FREEZESSRIGHT    => "FreezeStepSizeRightBound",
  :MASSMATRIX       => "MassMatrix",
  :JACOBIMATRIX     => "JacobiMatrix",
  :JACOBIBANDSTRUCT => "JacobiBandStructure",

  :MAXSTAGES        => "MaximalNumberOfStages",
  :MINSTAGES        => "MinimalNumberOfStages",
  :INITSTAGES       => "InitialNumberOfStages",
  :ORDERINCFACTOR   => "OrderIncreaseFactor",
  :ORDERDECFACTOR   => "OrderDecreaseFactor",
  :ORDERDECSTEPFAC1 => "OrderDecreaseStepFactor1",
  :ORDERDECSTEPFAC2 => "OrderDecreaseStepFactor2",

  :RHSAUTONOMOUS    => "AutonomousRHS",
  :LAMBDADENSE      => "LambdaForDenseOutput",
  :WORKFORRHS       => "WorkForRightHandSide",
  :WORKFORJAC       => "WorkForJacobimatrix",
  :WORKFORDEC       => "WorkForLuDecomposition",
  :WORKFORSOL       => "WorkForSubstitution",

  :BVPCLASS         => "BoundaryValueProblemClass",
  :SOLMETHOD        => "SolutionMethod",
  :IVPOPT           => "OptionsForIVPsolver")

function buildOptions(o,optionlist,aliases,aliases_reversed)
  dict1 = Dict{Symbol,Any}([Pair(k,o[k]) for k in (keys(o) ∩ optionlist)])
  dict2 = Dict([Pair(aliases_reversed[k],o[k]) for k in (keys(o) ∩ values(aliases))])
  merge(dict1,dict2)
end

function determine_initΔt(u₀,abstol,reltol,internalNorm,f,order)
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
