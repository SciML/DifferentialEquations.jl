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
              Δtmax=nothing,Δtmin=nothing,tType=typeof(Δt),internalNorm = 2,
              progressBar=false,progressSteps=1000)
  tspan = vec(tspan)
  if tspan[2]-tspan[1]<0 || length(tspan)>2
    error("tspan must be two numbers and final time must be greater than starting time. Aborting.")
  end

  @unpack prob:   f,u₀,knownSol,sol, numVars, sizeu

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

  if alg ∈ DIFFERENTIALEQUATIONSJL_ALGORITHMS
    if alg==:Euler
      u,t,uFull,tFull = ode_euler(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
    elseif alg==:Midpoint
      u,t,uFull,tFull = ode_midpoint(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
    elseif alg==:RK4
      u,t,uFull,tFull = ode_rk4(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,progressBar)
    elseif alg==:ExplicitRK
      u,t,uFull,tFull = ode_explicitrk(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,A,c,α,αEEst,stages,order,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar)
    elseif alg==:ImplicitEuler
      u,t,uFull,tFull = ode_impliciteuler(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu,progressBar)
    elseif alg==:Trapezoid
      u,t,uFull,tFull = ode_trapezoid(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu,progressBar)
    elseif alg==:Rosenbrock32
      u,t,uFull,tFull = ode_rosenbrock32(f,u,t,Δt,T,iter,maxIters,uFull,tFull,saveSteps,fullSave,adaptive,sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalNorm,progressBar)
    end
  elseif alg ∈ ODEINTERFACE_ALGORITHMS
    if typeof(u) <: Number
      u = [u]
    end
    @eval begin
      import ODEInterface
      export ODEInterface
      ODEInterface.loadODESolvers()
    end
    function outputfcn(reason,told,t,x,eval_sol_func,extra_data)
      #(progressBar && atomLoaded && iter%progressSteps==0) ? Main.Atom.progress(t/T) : nothing
      #println("outputfcn called with reason=$reason, told=$told, t=$t, x=$x")
      return ODEInterface.OUTPUTFCN_RET_CONTINUE
    end
    opt=ODEInterface.OptionsODE(
      ODEInterface.OPT_RTOL       => reltol,
      ODEInterface.OPT_ATOL       => abstol,
      ODEInterface.OPT_OUTPUTFCN  => outputfcn,
      ODEInterface.OPT_OUTPUTMODE => ODEInterface.OUTPUTFCN_WODENSE,
      ODEInterface.OPT_MAXSTEPS   => maxIters,
      ODEInterface.OPT_MAXSS      => Δtmax
      )
    if alg==:dopri5
      tFull,uFull,retcode,stats = ODEInterface.odecall(ODEInterface.dopri5,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opt)
    elseif alg==:dopri853
      tFull,uFull,retcode,stats = ODEInterface.odecall(ODEInterface.dop853,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opt)
    elseif alg==:odex
      tFull,uFull,retcode,stats = ODEInterface.odecall(ODEInterface.odex,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opt)
    elseif alg==:seulex
      tFull,uFull,retcode,stats = ODEInterface.odecall(ODEInterface.seulex,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opt)
    elseif alg==:radau
      tFull,uFull,retcode,stats = ODEInterface.odecall(ODEInterface.radau,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opt)
    elseif alg==:radau5
      tFull,uFull,retcode,stats = ODEInterface.odecall(ODEInterface.radau5,(t,u)->vec(f(reshape(u,sizeu),t)),Float64[t,T],vec(u),opt)
    end
    t = tFull[end]
    uFull = reshape(uFull,(length(tFull),sizeu...))
    if typeof(u₀)<:AbstractArray
      u[:] = vec(uFull[end,:,:]) #Quick hack
    else
      u = uFull[end]
    end
  elseif alg ∈ ODEJL_ALGORITHMS
    if typeof(u) <: Number
      u = [u]
      uFull = GrowableArray(u)
    end
    @eval begin
      import ODE
      export ODE
    end
    ode  = ODE.ExplicitODE(t,u,(t,y,dy)->dy[:]=f(y,t)) #not really inplace
    opts = Dict(:initstep=>Δt,
            :tstop=>T,
            #:tout=>[0.,0.5,1.],
            #:points=>:specified,
            :reltol=>reltol,
            :abstol=>abstol)

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
    out = ODE.solve(ode,stepper;opts...)
    for (t,u) in out # iterate over the solution
      push!(tFull,t)
      push!(uFull,u)
    end
    t = tFull[end]
    u = uFull[end]
    uFull = copy(uFull)
  end

  println(tFull)
  println(uFull)
  println(u₀)
  println(t)
  if knownSol
    uTrue = sol(u₀,t)
    if fullSave
      solFull = GrowableArray(sol(u₀,tFull[1]))
      for i in 2:size(uFull,1)
        push!(solFull,sol(u₀,tFull[i]))
      end
      uFull = copy(uFull)
      solFull = copy(solFull)
      #println(uFull)
      #println(solFull)
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

const DIFFERENTIALEQUATIONSJL_ALGORITHMS = Set([:Euler,:Midpoint,:RK4,:ExplicitRK,:ImplicitEuler,:Trapezoid,:Rosenbrock32])
const ODEINTERFACE_ALGORITHMS = Set([:dopri5,:dopri853,:odex,:radau5,:radau,:seulex])
const ODEJL_ALGORITHMS = Set([:ode23,:ode45,:ode78,:ode23s,:ode1,:ode2_midpoint,:ode2_heun,:ode4,:ode45_fe])
