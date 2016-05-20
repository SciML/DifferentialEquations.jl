function solve(sdeProb::SDEProblem,Δt,T;fullSave::Bool = false,saveSteps::Int = 1,alg::AbstractString="EM")
  @unpack sdeProb: f,σ,u₀,knownSol,sol
  N = ceil(Int64,T/Δt)

  u = float(u₀)
  t = 0.0
  W = 0.0
  if fullSave
    uFull = Vector{Float64}(ceil(Int64,N/saveSteps)+1)
    tFull = Vector{Float64}(ceil(Int64,N/saveSteps)+1)
    WFull = Vector{Float64}(ceil(Int64,N/saveSteps)+1)
    saveIdx = 1
    uFull[saveIdx] = u
    tFull[saveIdx] = t
    WFull[saveIdx] = W
  end

  sqΔt = sqrt(Δt)

  for i = 1:N
    dW = sqΔt*randn()
    u = u + Δt*f(u,t) + σ(u,t)*dW
    t=t+Δt
    W = W + dW
    if fullSave && i%saveSteps==0
      saveIdx+=1
      uFull[saveIdx] = u
      tFull[saveIdx] = t
      WFull[saveIdx] = W
    end
  end

  if knownSol
    uTrue = sol(u₀,t,W)
    if fullSave
      solFull = sol(u₀,tFull,WFull)
      return(SDESolution(u,uTrue,uFull=uFull,tFull=tFull,WFull=WFull,solFull=solFull))
    else
      return(SDESolution(u,uTrue))
    end
  else #No known sol
    return(SDESolution(u))
  end
end
