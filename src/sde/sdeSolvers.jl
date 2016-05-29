@def SRI begin
  chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
  chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
  chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h

  H0[:]=zeros(size(u)...,length(α))
  H1[:]=zeros(size(u)...,length(α))
  for i = 1:length(α)
    if numVars == 1
      A0temp = 0.0
      B0temp = 0.0
      A1temp = 0.0
      B1temp = 0.0
    else
      A0temp = zeros(size(u))
      B0temp = zeros(size(u))
      A1temp = zeros(size(u))
      B1temp = zeros(size(u))
    end
    for j = 1:i-1
      @inbounds A0temp += A₀[i,j]*f(H0[..,j],t + c₀[j]*Δt)
      @inbounds B0temp += B₀[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
      @inbounds A1temp += A₁[i,j]*f(H0[..,j],t + c₀[j]*Δt)
      @inbounds B1temp += B₁[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
    end
    H0[..,i] = u + A0temp*Δt + B0temp.*chi2
    H1[..,i] = u + A1temp*Δt + B1temp*sqΔt
  end
  if numVars == 1
    atemp = 0.0
    btemp = 0.0
  else
    atemp = zeros(size(u))
    btemp = zeros(size(u))
  end
  for i = 1:length(α)
    @inbounds atemp += α[i]*f(H0[..,i],t+c₀[i]*Δt)
    @inbounds btemp += (β₁[i]*ΔW + β₂[i]*chi1 + β₃[i]*chi2 + β₄[i]*chi3).*σ(H1[..,i],t+c₁[i]*Δt)
  end
  u = u + Δt*atemp + btemp
end

@def SRIVectorized begin
  chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
  chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
  chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h
  H0[:]=zeros(length(α))
  H1[:]=zeros(length(α))
  for i = 1:length(α)
    H0temp = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H1,t+c₁*Δt))
    H1[i]  = u + Δt*dot(vec(A₁[i,:]),f(H0,t + c₀*Δt)) + sqΔt*dot(vec(B₁[i,:]),σ(H1,t+c₁*Δt))
    H0[i] = H0temp
  end
  u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW + β₂*chi1 + β₃*chi2 + β₄*chi3,σ(H1,t+c₁*Δt))
end

@def SRAVectorized begin
  chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
  H0[:]=zeros(length(α))
  for i = 1:length(α)
    H0[i] = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H0,t+c₁*Δt))
  end
  u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW + β₂*chi2,σ(H0,t+c₁*Δt))
end

function solve(sdeProb::SDEProblem,Δt,T;fullSave::Bool = false,saveSteps::Int = 1,alg::AbstractString="EM")

  @unpack sdeProb: f,σ,u₀,knownSol,sol, numVars, sizeu

  u = float(u₀)
  t = 0.0
  if numVars == 1
    W = 0.0
    Z = 0.0
  else
    W = zeros(sizeu)
    Z = zeros(sizeu)
  end

  if fullSave
    uFull = GrowableArray(u)
    tFull = Vector{Float64}(0)
    WFull = GrowableArray(W)
    push!(tFull,t)
  end

  #PreProcess
  sqΔt = sqrt(Δt)

  if alg=="SRI"
    SRI = constructSRIW1()
    @unpack SRI: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
    if numVars == 1
      H0 = Array{Float64}(length(α))
      H1 = Array{Float64}(length(α))
    else
      H0 = Array{Float64}(size(u)...,length(α))
      H1 = Array{Float64}(size(u)...,length(α))
    end
  elseif alg=="SRA"
    SRA = constructSRA1()
    @unpack SRA: c₀,c₁,A₀,B₀,α,β₁,β₂
    if numVars == 1
      H0 = Array{Float64}(length(α))
    else
      H0 = Array{Float64}(size(u)...,length(α))
    end
  end
  if numVars == 1
    rands = ChunkedArray(randn)
  else
    rands = ChunkedArray(randn,u)
  end
  iter = 0
  while t < T
    iter += 1
    ΔW = sqΔt*next(rands)
    ΔZ = sqΔt*next(rands)
    if alg=="EM"
      u = u + Δt.*f(u,t) + σ(u,t).*ΔW
    elseif alg=="RKMil"
      K = u + Δt.*f(u,t)
      L = σ(u,t)
      utilde = K + L.*sqΔt
      u = K+L.*ΔW+(σ(utilde,t)-σ(u,t))./(2sqΔt).*(ΔW.^2 - Δt)
    elseif alg=="SRA" && numVars == 1
      @SRAVectorized
    elseif alg=="SRI" && numVars > 1 #Only for explicit
      @SRI
    elseif alg=="SRI" && numVars == 1 #Only for explicit
      @SRIVectorized
    end
    t = t + Δt
    W = W + ΔW
    Z = Z + ΔZ
    if fullSave && iter%saveSteps==0
      push!(uFull,u)
      push!(tFull,t)
      push!(WFull,W)
    end
  end

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
      return(SDESolution(u,uTrue,uFull=uFull,tFull=tFull,WFull=WFull,solFull=solFull))
    else
      return(SDESolution(u,uTrue))
    end
  else #No known sol
    return(SDESolution(u))
  end
end
