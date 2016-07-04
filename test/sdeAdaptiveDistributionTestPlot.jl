addprocs(CPU_CORES)

@everywhere begin
  using DifferentialEquations, Stats, Distributions, HypothesisTests,
        EllipsisNotation, Plots, JLD
  prob = linearSDEExample(α=1/10,β=1/20,u₀=1/2)
  tols = [1e-1;1e-3;1e-5]
  srand(199 + myid())
  T = 2
  N = 200
  M = 20
  K = length(tols)
  ps = Array{Float64}(M,K,3)
end

for k in eachindex(tols)
  println("k = $k")
  tol = tols[k]
  DifferentialEquations.sendto(workers(),tol=tol)
  println("RSwM1")
  @progress for j = 1:M
    println("j = $j")
    Wtmp = pmap((i)->(solve(prob::SDEProblem,[0,T],Δt=1/2^(4),alg="SRIW1Optimized",adaptive=true,abstol=tol,reltol=0,adaptiveAlg="RSwM1").W),1:N)
    Wends = Float64[float(Wtmp[i]) for i in 1:N]
    kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
    ps[j,k,1] = pvalue(kssol) #Should be not significant (most of the time)
  end
  println("The number rejected is $(sum(ps[:,k,1] .< 0.05))")
end

for k in eachindex(tols)
  println("k = $k")
  tol = tols[k]
  DifferentialEquations.sendto(workers(),tol=tol)
  println("RSwM2")
  @progress for j = 1:M
    println("j = $j")
    Wtmp = pmap((i)->solve(prob::SDEProblem,[0,T],Δt=1/2^(4),alg="SRIW1Optimized",adaptive=true,abstol=tol,reltol=0,adaptiveAlg="RSwM2").W,1:N)
    Wends = Float64[Wtmp[i] for i in 1:N]
    kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
    ps[j,k,2] = pvalue(kssol) #Should be not significant (most of the time)
  end
  println("The number rejected is $(sum(ps[:,k,2] .< 0.05))")
end

for k in eachindex(tols)
  println("k = $k")
  tol = tols[k]
  DifferentialEquations.sendto(workers(),tol=tol)
  println("RSwM3")
  @progress for j = 1:M
    println("j = $j")
    Wtmp = pmap((i)->solve(prob::SDEProblem,[0,T],Δt=1/2^(4),alg="SRIW1Optimized",adaptive=true,abstol=tol,reltol=0,adaptiveAlg="RSwM3").W,1:N)
    Wends = Float64[Wtmp[i] for i in 1:N]
    kssol = ApproximateOneSampleKSTest(Wends/sqrt(T), Normal())
    ps[j,k,3] = pvalue(kssol) #Should be not significant (most of the time)
  end
  println("The number rejected is $(sum(ps[:,k,3] .< 0.05))")
end

titleSize = 28
guideSize = 26
tickSize = 26
circleSize = 12
p = Vector{Any}(3)

p[1] = Plots.plot(vec(repmat(tols,1,M)'),vec(ps[..,1]),xscale=:log,yscale=:log,
      linetype=:scatter,yguide="Kolmogorov Smirnov P-values",title="RSwM1",left_margin=150px,
      top_margin=50px,xticks=[1e-1,1e-3,1e-5],bottom_margin=50px,markersize=circleSize,
      guidefont=font(guideSize),titlefont=font(titleSize),tickfont=font(tickSize),ylim=(1e-2,1),xlim=(1e-6,1),leg=false)

p[2] = Plots.plot(vec(repmat(tols,1,M)'),vec(ps[..,2]),xscale=:log10,yscale=:log10,linetype=:scatter,
      ylim=(1e-2,1),xlim=(1e-6,1),title="RSwM2",xguide="Absolute Tolerance",
      top_margin=50px,xticks=[1e-1,1e-3,1e-5],bottom_margin=50px,markersize=circleSize,
      guidefont=font(guideSize),titlefont=font(titleSize),tickfont=font(tickSize),ylim=(1e-2,1),xlim=(1e-6,1),leg=false)

p[3] = Plots.plot(vec(repmat(tols,1,M)'),vec(ps[..,3]),xscale=:log10,yscale=:log10,
        linetype=:scatter,ylim=(1e-2,1),xlim=(1e-6,1),
        top_margin=50px,xticks=[1e-1,1e-3,1e-5],bottom_margin=50px,
      title="RSwM3",right_margin=30px,markersize=circleSize,
      guidefont=font(guideSize),titlefont=font(titleSize),tickfont=font(tickSize),ylim=(1e-2,1),xlim=(1e-6,1),leg=false)

Plots.plot(p[1],p[2],p[3],size=(1200,800),layout=(1,3),link=:y)
Plots.gui()
