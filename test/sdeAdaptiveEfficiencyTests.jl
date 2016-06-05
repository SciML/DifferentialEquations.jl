using DifferentialEquations, Plots, LaTeXStrings


## Setup
N = 5
msims = Vector{MonteCarloSimulation}(N)
elapsed = Vector{Float64}(N)
medians = Vector{Float64}(N)
means   = Vector{Float64}(N)
tols    = Vector{Float64}(N)

## Problem 1
prob = linearSDEExample(α=1/10,β=1/20,u₀=1/2)
## Problem 2
#prob = waveSDEExample()
## Problem 3
#prob = additiveSDEExample()

adaptiveAlg = "RSwM1"
@progress for i=0:N-1
  tols[i+1] = 10.0^(-i)
  msims[i+1] = monteCarloSim(prob::SDEProblem,Δt=1//2^(4),adaptive=true,numMonte=1000,abstol=10.0^(-i),reltol=0,adaptiveAlg=adaptiveAlg)
  elapsed[i+1] = msims[i+1].elapsedTime
  medians[i+1] = msims[i+1].medians["l2"]
  means[i+1]   = msims[i+1].means["l2"]
end

l2 = L"l^2"
p11 = plot(tols,means,xscale=:log10,yscale=:log10,title="Test",titlefont=font(24),legendfont=font(14),tickfont=font(16),guidefont=font(16),xguide=L"ϵ",yguide="$l2 Error")
p12 = plot(tols,medians,xscale=:log10,yscale=:log10)
p13 = plot(tols,elapsed,xscale=:log10,yscale=:log10)

adaptiveAlg = "RSwM2"
@progress for i=0:N-1
  tols[i+1] = 10.0^(-i)
  msims[i+1] = monteCarloSim(prob::SDEProblem,Δt=1//2^(4),adaptive=true,numMonte=1000,abstol=10.0^(-i),reltol=0,adaptiveAlg=adaptiveAlg)
  elapsed[i+1] = msims[i+1].elapsedTime
  medians[i+1] = msims[i+1].medians["l2"]
  means[i+1]   = msims[i+1].means["l2"]
end

p21 = plot(tols,means,xscale=:log10,yscale=:log10)
p22 = plot(tols,medians,xscale=:log10,yscale=:log10)
p23 = plot(tols,elapsed,xscale=:log10,yscale=:log10)

adaptiveAlg = "RSwM3"
@progress for i=0:N-1
  tols[i+1] = 10.0^(-i)
  msims[i+1] = monteCarloSim(prob::SDEProblem,Δt=1//2^(4),adaptive=true,numMonte=1000,abstol=10.0^(-i),reltol=0,adaptiveAlg=adaptiveAlg)
  elapsed[i+1] = msims[i+1].elapsedTime
  medians[i+1] = msims[i+1].medians["l2"]
  means[i+1]   = msims[i+1].means["l2"]
end

p31 = plot(tols,means,xscale=:log10,yscale=:log10)
p32 = plot(tols,medians,xscale=:log10,yscale=:log10)
p33 = plot(tols,elapsed,xscale=:log10,yscale=:log10)

plot(p11,p21,p31,title="Means vs Tols")
plot(p12,p22,p32,title="Medians vs Tols")
plot(p13,p23,p33,title="Elapsed Time vs Tols")
