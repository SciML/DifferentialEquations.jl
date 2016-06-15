addprocs(CPU_CORES)

@everywhere begin
  using DifferentialEquations, Plots, EllipsisNotation, JLD, LaTeXStrings
  srand(99 + myid())
  prob = oval2ModelExample(largeFluctuations=true,useBigs=false)
  sol = solve(prob::SDEProblem,[0;1],Δt=1/2^(8),fullSave=true,alg="EM",adaptive=false,progressBar=true,saveSteps=1,abstol=1e-6,reltol=1e-4)
  Int(sol.u[1]!=NaN)
  tspan = [0;1]
  js = 18:18
  Δts = 1./2.^(js)
  fails = Array{Int}(length(Δts),2)
  times = Array{Float64}(length(Δts),2)
  numRuns = 100
end

for j in eachindex(js)
  println("j = $j")
  function runEM(i,j)
    sol =solve(prob::SDEProblem,tspan,Δt=Δts[j],alg="EM")
    Int(any(isnan,sol.u))
  end
  t1 = @elapsed numFails = sum(pmap((i)->runEM(i,j),1:numRuns))
  println("The number of Euler-Maruyama Fails is $numFails. Elapsed time was $t1")
  fails[j,1] = numFails
  times[j,1] = t1
end


for j in eachindex(js)
  println("j = $j")
  function runSRI(i,j)
    sol =solve(prob::SDEProblem,tspan,Δt=Δts[j],alg="SRIW1Optimized")
    Int(any(isnan,sol.u))
  end
  t2 = @elapsed numFails = sum(pmap((i)->runSRI(i,j),1:numRuns))
  println("The number of Rossler-SRI Fails is $numFails. Elapsed time was $t2")
  fails[j,2] = numFails
  times[j,2] = t2
end


@everywhere function runAdaptive(i)
  sol = solve(prob::SDEProblem,tspan,Δt=1/2^(8),alg="SRIW1Optimized",adaptiveAlg="RSwM3",adaptive=true,abstol=1e-5,reltol=1e-3)
  Int(any(isnan,sol.u))
end
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptive,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")



lw = 3
p2 = plot(Δts,times,xscale=:log2,yscale=:log2,guidefont=font(16),tickfont=font(14),yguide="Elapsed Time (s)",xguide=L"Chosen $\Delta t$",top_margin=50px,linewidth=lw,lab=["Euler-Maruyama" "SRIW1"],legendfont=font(14))
plot!(Δts,repmat([adaptiveTime],11),linewidth=lw,line=:dash,lab="ESRK+RSwM3",left_margin=75px)
plot(p2,size=(800,800))
gui()




#Big Run
@time solve(prob::SDEProblem,[0;10],Δt=(1/2)^(8),fullSave=true,alg="SRIW1Optimized",
          adaptiveAlg="RSwM3",adaptive=true,progressBar=true,
          saveSteps=100,abstol=1e-5,reltol=1e-3)

#Plots
lw = 2

p1 = plot(sol.tFull,sol.uFull[..,16],top_margin=50px,title="Ecad",xguide="Time",yguide="Concentration",guidefont=font(16),tickfont=font(16),linewidth=lw,left_margin=85px,leg=false)
p2 = plot(sol.tFull,sol.uFull[..,17],top_margin=50px,title="Vim",xguide="Time",yguide="Concentration",guidefont=font(16),tickfont=font(16),linewidth=lw,leg=false)
p3 = plot(sol.tFull,sol.ΔtFull,xguide="Time",yguide="Accepted Dt",guidefont=font(16),tickfont=font(16),yscale=:log10,linewidth=lw,left_margin=110px,bottom_margin=65px,leg=false)
plot(p1,p2,p3,layout=@layout([a b;c]),size=(1200,800),title="Adaptive Solution to Stochastic Cell Model")
gui()

p1 = plot(sol.tFull,sol.uFull[..,16],top_margin=50px,title="Ecad",xguide="Time",yguide="Concentration",guidefont=font(16),tickfont=font(16))
p2 = plot(sol.tFull,sol.uFull[..,17],top_margin=50px,title="Vim",xguide="Time",yguide="Concentration",guidefont=font(16),tickfont=font(16))
p3 = plot(sol.tFull[1:end-1],diff(sol.tFull)/100,xguide="Time",yguide="Accepted Dt 100 Step Averages",guidefont=font(16),tickfont=font(16),yscale=:log10)
plot(p1,p2,p3,layout=@layout([a b;c]),size=(1200,800),title="Adaptive Solution to Stochastic Cell Model")
gui()

save("Oval2Solution.jld","sol",sol,"prob",prob)
