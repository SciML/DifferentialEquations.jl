using DifferentialEquations, Stats, Distributions, HypothesisTests,
      EllipsisNotation, Plots, JLD
prob = linearSDEExample()
tols = [1e-1;1e-3;1e-5]
srand(200)
T = 1
N = 200
M = 20
K = length(tols)
ps = Array{Float64}(M,K,3)

for k in eachindex(tols)
  println("k = $k")
  tol = tols[k]
  println("RSwM1")
  @progress for j = 1:M
    Wends = Vector{Float64}(N)
    for i = 1:N
      sol =solve(prob::SDEProblem,[0,T],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=tol,reltol=0,adaptiveAlg="RSwM1")
      Wends[i] = sol.WFull[end]
    end
    kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
    ps[j,k,1] = pvalue(kssol) #Should be not significant (most of the time)
  end

  println("RSwM2")
  @progress for j = 1:M
    Wends = Vector{Float64}(N)
    for i = 1:N
      sol =solve(prob::SDEProblem,[0,T],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=tol,reltol=0,adaptiveAlg="RSwM2")
      Wends[i] = sol.WFull[end]
    end
    kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
    ps[j,k,2] = pvalue(kssol) #Should be not significant (most of the time)
  end

  println("RSwM3")
  @progress for j = 1:M
    Wends = Vector{Float64}(N)
    for i = 1:N
      sol =solve(prob::SDEProblem,[0,T],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=tol,reltol=0,adaptiveAlg="RSwM3")
      Wends[i] = sol.WFull[end]
    end
    kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
    ps[j,k,3] = pvalue(kssol) #Should be not significant (most of the time)
  end
end

tol = 1e-4
p = Vector{Float64}(M)
@progress for j = 1:M
  Wends = Vector{Float64}(N)
  for i = 1:N
    sol =solve(prob::SDEProblem,[0,T],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=tol,reltol=0,adaptiveAlg="RSwM3")
    Wends[i] = sol.WFull[end]
  end
  kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
  p[j] = pvalue(kssol)
end
pvalue(kssol)>0.05

p = Vector{Any}(3)

p[1] = Plots.plot(vec(repmat(tols,1,M)'),vec(ps[..,1]'),xscale=:log10,yscale=:log10,linetype=:scatter,yguide="Kolmogorov Smirnov P-values",ylim=(1e-3,2),xlim=(1e-4,0),title="RSwM1",left_margin=90px,guidefont=font(16),titlefont=font(20),tickfont=font(16))

p[2] = Plots.plot(vec(repmat(tols,1,M)'),vec(ps[..,2]'),xscale=:log10,yscale=:log10,linetype=:scatter,ylim=(1e-3,2),xlim=(1e-4,0),title="RSwM2",xguide="Absolute Tolerance",guidefont=font(16),titlefont=font(20),tickfont=font(16))

p[3] = Plots.plot(vec(repmat(tols,1,M)'),vec(ps[..,3]'),xscale=:log10,yscale=:log10,linetype=:scatter,
        ylim=(1e-3,2),xlim=(1e-4,0),title="RSwM3",right_margin=30px,titlefont=font(20),tickfont=font(16))

Plots.plot(p[1],p[2],p[3],size=(1200,800),layout=(1,3))
Plots.gui()

save("ksTestResults.jld","ps",ps)
