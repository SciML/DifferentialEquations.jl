function solve(prob::ODEProblem,Δt,T;fullSave::Bool = false,saveSteps::Int = 1,alg::AbstractString="Euler",tableau=DEFAULT_TABLEAU,adaptive=false,γ=2,tol=1e-4,qmax=10)

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
        absEEst = norm(utilde-uEEst,2)
        relEEst = absEEst/norm(uEEst,2)
      else
        u = u + Δt*utilde
      end
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
