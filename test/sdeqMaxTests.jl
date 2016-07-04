addprocs(CPU_CORES)

@everywhere qs = 1.0 + 2.0.^(-5:2)
times = Array{Float64}(length(qs),4)
means = Array{Float64}(length(qs),4)

@everywhere begin
  using DifferentialEquations, Plots, EllipsisNotation, LaTeXStrings
  srand(99 + myid())
  prob = oval2ModelExample(largeFluctuations=true,useBigs=false)
  println("Solve once to compile.")
  sol = solve(prob::SDEProblem,[0;1],Δt=1/2^(18),fullSave=true,alg="EM",adaptive=false)
  Int(sol.u[1]!=NaN)
  println("Compilation complete.")
  tspan = [0;1]
  numRuns = 10000

  probs = Vector{SDEProblem}(3)
  p1 = Vector{Any}(3)
  p2 = Vector{Any}(3)
  p3 = Vector{Any}(3)
  ## Problem 1
  probs[1] = linearSDEExample(α=1/10,β=1/20,u₀=1/2)
  ## Problem 2
  probs[2] = waveSDEExample()
  ## Problem 3
  probs[3] = additiveSDEExample()
end
println("Setup Complete")






## Timing Runs

@everywhere function runAdaptive(i,k)
  sol = solve(prob::SDEProblem,tspan,Δt=1/2^(8),alg="SRIW1Optimized",adaptiveAlg="RSwM3",adaptive=true,abstol=2.0^(-14),reltol=2.0^(-9),maxIters=Int(1e12),qmax=qs[k])
  Int(any(isnan,sol.u))
end


for k in eachindex(qs)
  DifferentialEquations.sendto(workers(), k=k)
  @everywhere srand(99 + myid())
  adaptiveTime = @elapsed numFails = sum(pmap((i)->runAdaptive(i,k),1:numRuns))
  println("k was $k. The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")
  times[k,4] = adaptiveTime
end

####################### Other Examples
#Compile
monteCarloSim(probs[1]::SDEProblem,Δt=1/2^(4),adaptive=true,numMonte=1000,abstol=2.0^(-1),reltol=0,adaptiveAlg="RSwM3",alg="SRIW1Optimized")

@everywhere numRuns*=10
for k in eachindex(probs)
  println("Problem $k")
  ## Setup
  prob = probs[k]
  DifferentialEquations.sendto(workers(), prob=prob)

  for i in eachindex(qs)
    DifferentialEquations.sendto(workers(), i=i)
    msim = monteCarloSim(prob::SDEProblem,Δt=1/2^(4),adaptive=true,numMonte=numRuns,abstol=2.0^(-13),reltol=0,adaptiveAlg="RSwM3",alg="SRIW1Optimized",qmax=qs[i])
    times[i,k] = msim.elapsedTime
    means[i,k] = msim.means["final"]
    println("for k=$k and i=$i, we get that the error was $(means[i,k]) and it took $(times[i,k]) seconds")
  end
end
