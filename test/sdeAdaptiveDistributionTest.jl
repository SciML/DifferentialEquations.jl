using DifferentialEquations, Stats, Distributions, HypothesisTests
prob = linearSDEExample()
srand(200)
T = 1
N = 100
M=5
ps = Vector{Float64}(M)


for j = 1:M
  Wends = Vector{Float64}(N)
  for i = 1:N
    sol =solve(prob::SDEProblem,[0,T],Δt=1/2^(4),fullSave=true,alg=:SRI,adaptive=true,abstol=1e-2,reltol=0,adaptivealg=:RSwM1)
    Wends[i] = sol.WFull[end]
  end
  kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
  ps[j] = pvalue(kssol) #Should be not significant (most of the time)
end

bool1 = sum(ps .> 0.05) > length(ps)/2 ### Make sure more passes than fails

for j = 1:M
  Wends = Vector{Float64}(N)
  for i = 1:N
    sol =solve(prob::SDEProblem,[0,T],Δt=1/2^(4),fullSave=true,alg=:SRI,adaptive=true,abstol=1e-2,reltol=0,adaptivealg=:RSwM2)
    Wends[i] = sol.WFull[end]
  end
  kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
  ps[j] = pvalue(kssol) #Should be not significant (most of the time)
end

bool2 = sum(ps .> 0.05) > length(ps)/2 ### Make sure more passes than fails

for j = 1:M
  Wends = Vector{Float64}(N)
  for i = 1:N
    sol =solve(prob::SDEProblem,[0,T],Δt=1/2^(4),fullSave=true,alg=:SRI,adaptive=true,abstol=1e-2,reltol=0,adaptivealg=:RSwM3)
    Wends[i] = sol.WFull[end]
  end
  kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
  ps[j] = pvalue(kssol) #Should be not significant (most of the time)
end

bool3 = sum(ps .> 0.05) > length(ps)/2 ### Make sure more passes than fails
