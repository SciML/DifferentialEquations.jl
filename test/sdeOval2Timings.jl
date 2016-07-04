addprocs(CPU_CORES)

@everywhere begin
  using DifferentialEquations, Plots, EllipsisNotation, LaTeXStrings
  srand(99 + myid())
  prob = oval2ModelExample(largeFluctuations=true,useBigs=false)
  println("Solve once to compile.")
  sol = solve(prob::SDEProblem,[0;1],Δt=1/2^(18),fullSave=true,alg="EM",adaptive=false)
  Int(sol.u[1]!=NaN)
  println("Compilation complete.")
  tspan = [0;1]
  js = 16:21
  Δts = 1./2.^(js)
  fails = Array{Int}(length(Δts),3)
  times = Array{Float64}(length(Δts),3)
  numRuns = 10000
end
println("Setup Complete")

## Timing Runs

@everywhere function runAdaptive(i)
  sol = solve(prob::SDEProblem,tspan,Δt=1/2^(8),alg="SRIW1Optimized",adaptiveAlg="RSwM3",adaptive=true,abstol=2.0^(-14),reltol=2.0^(-9),maxIters=Int(1e12),qmax=1.125)
  Int(any(isnan,sol.u))
end
@everywhere srand(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptive,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


for j in eachindex(js)
  println("j = $j")
  DifferentialEquations.sendto(workers(), j=j)
  @everywhere function runEM(i,j)
    sol =solve(prob::SDEProblem,tspan,Δt=Δts[j],alg="EM",maxIters=Int(1e11))
    Int(any(isnan,sol.u))
  end
  @everywhere srand(99 + myid())
  t1 = @elapsed numFails = sum(pmap((i)->runEM(i,j),1:numRuns))
  println("The number of Euler-Maruyama Fails is $numFails. Elapsed time was $t1")
  fails[j,1] = numFails
  times[j,1] = t1
end


for j in 4:6
  println("j = $j")
  DifferentialEquations.sendto(workers(), j=j)
  @everywhere function runSRI(i,j)
    sol =solve(prob::SDEProblem,tspan,Δt=Δts[j],alg="SRIW1Optimized")
    Int(any(isnan,sol.u))
  end
  @everywhere srand(99 + myid())
  t2 = @elapsed numFails = sum(pmap((i)->runSRI(i,j),1:numRuns))
  println("The number of Rossler-SRI Fails is $numFails. Elapsed time was $t2")
  fails[j,2] = numFails
  times[j,2] = t2
end


for j in eachindex(js)
  println("j = $j")
  DifferentialEquations.sendto(workers(), j=j)
  @everywhere function runMil(i,j)
    sol =solve(prob::SDEProblem,tspan,Δt=Δts[j],alg="RKMil")
    Int(any(isnan,sol.u))
  end
  @everywhere srand(99 + myid())
  t3 = @elapsed numFails = sum(pmap((i)->runMil(i,j),1:numRuns))
  println("The number of RK-Milstein Fails is $numFails. Elapsed time was $t3")
  fails[j,3] = numFails
  times[j,3] = t3
end

lw = 3
p2 = plot(Δts,times,xscale=:log2,yscale=:log2,guidefont=font(16),tickfont=font(14),yguide="Elapsed Time (s)",xguide=L"Chosen $\Delta t$",top_margin=50px,linewidth=lw,lab=["Euler-Maruyama" "RK-Mil" "RosslerSRI"],legendfont=font(14))
plot!(Δts,repmat([adaptiveTime],11),linewidth=lw,line=:dash,lab="ESRK+RSwM3",left_margin=75px)
scatter!([2.0^(-20);2.0^(-20);2.0^(-18)],[times[5,1];times[5,2];times[3,3]],markersize=20,c=:red,lab="")
plot(p2,size=(800,800))
gui()
