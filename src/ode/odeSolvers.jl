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
* tableau - Takes in an object which defines a tableau. Default is Dormand-Prince 4/5.
* adaptive - Turns on adaptive timestepping for appropriate methods. Default is false.
* tol - The error tolerance of the adaptive method. Default is 1e-4.
* γ - The risk-factor γ in the q equation for adaptive timestepping. Default is 2.
* qmax - Defines the maximum value possible for the adaptive q. Default is 10.
"""
function solve(prob::ODEProblem,Δt::Number,T::Number;fullSave::Bool = false,saveSteps::Int = 1,alg::AbstractString="RK4",tableau=DEFAULT_TABLEAU,adaptive=false,γ=2,tol=1e-4,qmax=10)

  @unpack prob: f,u₀,knownSol,sol, numVars, sizeu
  u = float(u₀)
  t = 0.0

  if fullSave
    uFull = GrowableArray(u)
    tFull = Vector{Float64}(0)
    push!(tFull,t)
  end

  iter = 0

  #Pre-process
  if alg == "Midpoint"
    utilde = similar(u)
    halfΔt = .5Δt
  elseif alg == "RK4"
    k₁ = similar(u)
    k₂ = similar(u)
    k₃ = similar(u)
    k₄ = similar(u)
    halfΔt = .5Δt
  elseif alg == "ExplicitRK"
    # tableau from keyword argument
    @unpack tableau:   A,c,α,αEEst,stages,order
    ks = Array{Float64}(size(u)...,stages)
  elseif alg == "ImplicitEuler"
    function rhs(u,resid,uOld,t,Δt)
      u = reshape(u,sizeu...)
      resid = reshape(resid,sizeu...)
      resid[:] = u - uOld - Δt*f(u,t+Δt)
      vec(u)
      vec(resid)
    end
  elseif alg == "Trapezoid"
    function rhs(u,resid,uOld,t,Δt)
      u = reshape(u,sizeu...)
      resid = reshape(resid,sizeu...)
      resid[:] = u - uOld - .5Δt*(f(u,t+Δt)+f(uOld,t))
      u = vec(u)
      resid = vec(resid)
    end
  end

  while t < T
    iter += 1
    if alg=="Euler"
      u = u + Δt.*f(u,t)
    elseif alg=="Midpoint"
      utilde[:] = u + Δt.*f(u,t)
      u = u + Δt.*f(u+halfΔt*utilde,t+halfΔt)
    elseif alg=="RK4"
      k₁[:] = f(u,t)
      ttmp = t+halfΔt
      k₂[:] = f(u+halfΔt*k₁,ttmp)
      k₃[:] = f(u+halfΔt*k₂,ttmp)
      k₄[:] = f(u+Δt*k₃,t+Δt)
      u = u + Δt*(k₁ + 2k₂ + 2k₃ + k₄)/6
    elseif alg=="ExplicitRK"
      for i = 1:stages
        if typeof(u)<:Number
          utilde = 0
        else
          utilde = zeros(u)
        end
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
        absEEst = norm(utilde-uEEst,2)
        relEEst = absEEst/norm(uEEst,2)
      else
        u = u + Δt*utilde
      end
    elseif alg=="ImplicitEuler"
      uOld = copy(u)
      u = vec(u)
      nlres = nlsolve((u,resid)->rhs(u,resid,uOld,t,Δt),u)
      u = reshape(nlres.zero,sizeu...)
    elseif alg=="Trapezoid"
      uOld = copy(u)
      u = vec(u)
      nlres = nlsolve((u,resid)->rhs(u,resid,uOld,t,Δt),u)
      u = reshape(nlres.zero,sizeu...)
    end
    if adaptive
      standard = abs((tol)/(γ*absEEst)).^(1/order)
      if isinf(standard)
          q = qmax
      else
         q = min(qmax,max(standard,eps()))
      end
      if q > 1
        t = t + Δt
        u = utmp
        if fullSave && iter%saveSteps==0
          push!(uFull,u)
          push!(tFull,t)
        end
      end
      Δt = q*Δt
    else # Non adaptive
      t = t + Δt
      if fullSave && iter%saveSteps==0
        push!(uFull,u)
        push!(tFull,t)
      end
    end
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
    return(ODESolution(u))
  end
end
