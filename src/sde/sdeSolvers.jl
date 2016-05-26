function solve(sdeProb::SDEProblem,Δt,T;fullSave::Bool = false,saveSteps::Int = 1,alg::AbstractString="EM")

  @unpack sdeProb: f,σ,u₀,knownSol,sol
  N = ceil(Int,T/Δt)

  u = float(u₀)
  t = 0.0
  W = 0.0
  Z = 0.0
  if fullSave
    uFull = Vector{Float64}(ceil(Int,N/saveSteps)+1)
    tFull = Vector{Float64}(ceil(Int,N/saveSteps)+1)
    WFull = Vector{Float64}(ceil(Int,N/saveSteps)+1)
    saveIdx = 1
    uFull[saveIdx] = u
    tFull[saveIdx] = t
    WFull[saveIdx] = W
  end

  #PreProcess
  sqΔt = sqrt(Δt)

  if alg=="SRI"
    SRI = constructSRIW1()
    @unpack SRI: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
    H0 = Vector{Float64}(length(α))
    H1 = Vector{Float64}(length(α))
  elseif alg=="SRA"
    SRA = constructSRA1()
    @unpack SRA: c₀,c₁,A₀,B₀,α,β₁,β₂
    H0 = Vector{Float64}(length(α))
  end

  for i = 1:N
    ΔW = sqΔt*randn()
    ΔZ = sqΔt*randn()
    if alg=="EM"
      u = u + Δt*f(u,t) + σ(u,t)*ΔW
    elseif alg=="RKMil"
      #utilde = u + Δt*f(u,t) + sqdt*σ(u,t)
      #u = u + Δt*f(u,t) + σ(u,t)*ΔW + (utilde-)
      K = u + Δt*f(u,t)
      L = σ(u,t)
      utilde = K + L*sqΔt
      u = K+L*ΔW+(σ(utilde,t)-σ(u,t))./(2sqΔt).*(ΔW.^2 - Δt)
    elseif alg=="SRA"
      chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
      H0[:]=zeros(length(α))
      for i = 1:length(α)
        H0[i] = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H0,t+c₁*Δt))
      end
      u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW + β₂*chi2,σ(H0,t+c₁*Δt))
    elseif alg=="SRI"
      chi1 = .5*(ΔW^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
      chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
      chi3 = 1/6 * (ΔW^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h
      H0[:]=zeros(length(α))
      H1[:]=zeros(length(α))
      for i = 1:length(α)
        H0temp = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H1,t+c₁*Δt))
        H1[i]  = u + Δt*dot(vec(A₁[i,:]),f(H0,t + c₀*Δt)) + sqΔt*dot(vec(B₁[i,:]),σ(H1,t+c₁*Δt))
        H0[i] = H0temp
      end
      u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW + β₂*chi1 + β₃*chi2 + β₄*chi3,σ(H1,t+c₁*Δt))
    end
    t = t + Δt
    W = W + ΔW
    Z = Z + ΔZ
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
