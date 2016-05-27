function solve(sdeProb::SDEProblem,Δt,T;fullSave::Bool = false,saveSteps::Int = 1,alg::AbstractString="EM")

  @unpack sdeProb: f,σ,u₀,knownSol,sol, numVars, sizeu

  u = u₀
  t = 0.0
  W = zeros(sizeu)
  Z = zeros(sizeu)

  if fullSave
    uFull = Vector{typeof(u)}(0)
    tFull = Vector{Float64}(0)
    WFull = Vector{typeof(u)}(0)
    push!(uFull,u)
    push!(tFull,t)
    push!(WFull,W)
  end

  #PreProcess
  sqΔt = sqrt(Δt)

  if alg=="SRI"
    SRI = constructSRIW1()
    @unpack SRI: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
    H0 = Array{Float64}(length(α),numVars)
    H1 = Array{Float64}(length(α),numVars)
  elseif alg=="SRA"
    SRA = constructSRA1()
    @unpack SRA: c₀,c₁,A₀,B₀,α,β₁,β₂
    H0 = Array{Float64}(length(α),numVars)
  end

  i = 0
  while t < T
    i += 1
    ΔW = sqΔt*randn(sizeu)
    ΔZ = sqΔt*randn(sizeu)
    if alg=="EM"
      u = u + Δt.*f(u,t) + σ(u,t).*ΔW
    elseif alg=="RKMil"
      K = u + Δt.*f(u,t)
      L = σ(u,t)
      utilde = K + L.*sqΔt
      u = K+L.*ΔW+(σ(utilde,t)-σ(u,t))./(2sqΔt).*(ΔW.^2 - Δt)
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
      push!(uFull,u)
      push!(tFull,t)
      push!(WFull,W)
    end
  end

  if knownSol

    uTrue = sol(u₀,t,W)
    if fullSave
      ## Reshape solution. Needs to be fixed for arbitrary dimensions
      fill = Array{Float64}(size(WFull)...,size(u)...)
      fill2= Array{Float64}(size(WFull)...,size(u)...)
      for i in eachindex(uFull) ### Change to arbitrary dimensions
        fill[i,:] = WFull[i]'
        fill2[i,:]= uFull[i]'
      end
      uFull = fill2
      WFull = fill
      solFull = similar(WFull)
      for i in 1:size(WFull,1)
        solFull[i,:] = sol(u₀,tFull[i],WFull[i,:]')'
      end
      return(SDESolution(u,uTrue,uFull=uFull,tFull=tFull,WFull=WFull,solFull=solFull))
    else
      return(SDESolution(u,uTrue))
    end
  else #No known sol
    return(SDESolution(u))
  end
end
