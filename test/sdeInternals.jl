import Base.getindex, Base.setindex!
const .. = Val{:...}

setindex!{T}(A::AbstractArray{T,1}, x, ::Type{Val{:...}}, n) = A[n] = x
setindex!{T}(A::AbstractArray{T,2}, x, ::Type{Val{:...}}, n) = A[ :, n] = x
setindex!{T}(A::AbstractArray{T,3}, x, ::Type{Val{:...}}, n) = A[ :, :, n] =x

getindex{T}(A::AbstractArray{T,1}, ::Type{Val{:...}}, n) = A[n]
getindex{T}(A::AbstractArray{T,2}, ::Type{Val{:...}}, n) = A[ :, n]
getindex{T}(A::AbstractArray{T,3}, ::Type{Val{:...}}, n) = A[ :, :, n]

setindex!{T}(A::AbstractArray{T,1}, x, n, ::Type{Val{:...}}) = A[n] = x
setindex!{T}(A::AbstractArray{T,2}, x, n, ::Type{Val{:...}}) = A[n, :] = x
setindex!{T}(A::AbstractArray{T,3}, x, n, ::Type{Val{:...}}) = A[n, :, :] =x

getindex{T}(A::AbstractArray{T,1}, n, ::Type{Val{:...}}) = A[n]
getindex{T}(A::AbstractArray{T,2}, n, ::Type{Val{:...}}) = A[n, :]
getindex{T}(A::AbstractArray{T,3}, n, ::Type{Val{:...}}) = A[n, :, :]

using DifferentialEquations, Parameters
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
#sol =solve(prob::SDEProblem,1//2^(4),1,fullSave=true,alg="RKMil")

Δt = 1//2^4
T = 1
fullSave = false
saveSteps = 1
alg="SRI"

@unpack prob: f,σ,u₀,knownSol,sol, numVars, sizeu

u = u₀
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
  H0 = Array{Float64}(size(u)...,length(α))
  H1 = Array{Float64}(size(u)...,length(α))
  H02 = Array{Float64}(length(α))
  H12 = Array{Float64}(length(α))
elseif alg=="SRA"
  SRA = constructSRA1()
  @unpack SRA: c₀,c₁,A₀,B₀,α,β₁,β₂
  H0 = Array{Float64}(length(α),size(u)...)
end

iter = 0
while t < T
  iter += 1
  if numVars == 1
    ΔW = sqΔt*randn()
    ΔZ = sqΔt*randn()
  else
    ΔW = sqΔt*randn(sizeu)
    ΔZ = sqΔt*randn(sizeu)
  end
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
  elseif alg=="SRI" #Only for explicit
    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h
    H0[:]=zeros(size(u)...,length(α))
    H1[:]=zeros(size(u)...,length(α))
    for i = 1:length(α)
      A0temp = zeros(size(u))
      B0temp = zeros(size(u))
      A1temp = zeros(size(u))
      B1temp = zeros(size(u))
      for j = 1:i-1
        A0temp += A₀[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        B0temp += B₀[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
        A1temp += A₁[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        B1temp += B₁[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
      end
      H0[..,i] = u + A0temp*Δt + B0temp.*chi2
      H1[..,i] = u + A1temp*Δt + B1temp*sqΔt
    end
    #ab sum part is correct
    atemp = 0.0
    btemp = 0.0
    for i = 1:length(α)
      atemp += α[i]*f(H0[..,i],t+c₀[i]*Δt)
      btemp += (β₁[i]*ΔW + β₂[i]*chi1 + β₃[i]*chi2 + β₄[i]*chi3).*σ(H1[..,i],t+c₁[i]*Δt)
    end
    #=
    H02[:]=zeros(length(α))
    H12[:]=zeros(length(α))
    for i = 1:length(α)
      H02temp = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H1,t+c₁*Δt))
      H12[i]  = u + Δt*dot(vec(A₁[i,:]),f(H0,t + c₀*Δt)) + sqΔt*dot(vec(B₁[i,:]),σ(H1,t+c₁*Δt))
      H02[i] = H02temp
    end

    u2 = u + Δt*dot(α,f(H02,t+c₀*Δt)) + dot(β₁*ΔW + β₂*chi1 + β₃*chi2 + β₄*chi3,σ(H12,t+c₁*Δt))
    =#
    u = u + Δt*atemp + btemp
    #u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW + β₂*chi1 + β₃*chi2 + β₄*chi3,σ(H1,t+c₁*Δt))
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
