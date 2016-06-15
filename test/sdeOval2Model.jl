addprocs(32)
using DifferentialEquations, Plots, EllipsisNotation, JLD, LaTeXStrings
srand(100)
prob = oval2ModelExample(largeFluctuations=true,useBigs=false)




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

##Progress Bar

@everywhere function durationstring(nsec)
    days = div(nsec, 60*60*24)
    r = nsec - 60*60*24*days
    hours = div(r,60*60)
    r = r - 60*60*hours
    minutes = div(r, 60)
    seconds = r - 60*minutes

    hhmmss = @sprintf "%u:%02u:%02u" hours minutes seconds
    if days>0
       return @sprintf "%u days, %s" days hhmmss
    end
    hhmmss
end

##Adaptivity Necessity Tests
sol = solve(prob::SDEProblem,[0;1],Δt=1/2^(8),fullSave=true,alg="EM",adaptive=false,progressBar=true,saveSteps=1,abstol=1e-6,reltol=1e-4)
Int(sol.u[1]!=NaN)
tspan = [0;1]
js = 18:20
Δts = 1./2.^(js)
fails = Array{Int}(length(Δts),2)
times = Array{Float64}(length(Δts),2)
numRuns = 100

@progress for j in eachindex(js)
  println("j = $j")
  numFails = 0
  startTime = time()
  complete = SharedArray(Int,(numRuns),init=zeros(Int,numRuns))
  t1 = @elapsed numFails = @parallel (+) for i = 1:numRuns
    sol =solve(prob::SDEProblem,tspan,Δt=Δts[j],alg="EM")
    complete[i] = 1

    percentage_complete = 100*(sum(complete)/numRuns)
    elapsed_time = time() - startTime
    est_total_time = 100 * elapsed_time / percentage_complete
    eta_sec = round(Int, est_total_time - elapsed_time )
    eta = durationstring(eta_sec)
    println("\% Complete: $percentage_complete. ETA: $eta")
    flush(STDOUT)

    Int(any(isnan,sol.u))
  end
  fails[j,1] = numFails
  times[j,1] = t1
  println("The number of Euler-Maruyama Fails is $numFails. Elapsed time was $t1")
  numFails = 0
end


@progress for j in js
  numFails = 0
  complete = SharedArray(Int,(numRuns),init=zeros(Int,numRuns))
  startTime = time()
  t2 = @elapsed numFails = @parallel (+) for i = 1:numRuns
    sol =solve(prob::SDEProblem,tspan,Δt=Δts[j],alg="SRI",adaptive=false)
    complete[i] = 1

    percentage_complete = 100*(sum(complete)/numRuns)
    elapsed_time = time() - startTime
    est_total_time = 100 * elapsed_time / percentage_complete
    eta_sec = round(Int, est_total_time - elapsed_time )
    eta = durationstring(eta_sec)
    println("\% Complete: $percentage_complete. ETA: $eta")
    flush(STDOUT)
    Int(any(isnan,sol.u))
  end
  println("The number of Rossler-SRI Fails is $numFails. Elapsed time was $t2")
  fails[j,2] = numFails
  times[j,2] = t2
  numFails = 0
end

complete = SharedArray(Int,(numRuns),init=zeros(Int,numRuns))
startTime = time()
adaptiveTime = @elapsed numFails = @parallel (+) for i = 1:numRuns
  sol = solve(prob::SDEProblem,tspan,Δt=1/2^(8),alg="SRIW1Optimized",adaptiveAlg="RSwM3",adaptive=true,abstol=1e-5,reltol=1e-3)
  complete[i] = 1

  percentage_complete = 100*(sum(complete)/numRuns)
  elapsed_time = time() - startTime
  est_total_time = 100 * elapsed_time / percentage_complete
  eta_sec = round(Int, est_total_time - elapsed_time )
  eta = durationstring(eta_sec)
  println("\% Complete: $percentage_complete. ETA: $eta")
  flush(STDOUT)

  Int(any(isnan,sol.u))
end
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")



complete = SharedArray(Int,(numRuns),init=zeros(Int,numRuns))
startTime = time()
function runAdaptive(i)
  sol = solve(prob::SDEProblem,tspan,Δt=1/2^(8),alg="SRIW1Optimized",adaptiveAlg="RSwM3",adaptive=true,abstol=1e-5,reltol=1e-3)
  complete[i] = 1

  percentage_complete = 100*(sum(complete)/numRuns)
  elapsed_time = time() - startTime
  est_total_time = 100 * elapsed_time / percentage_complete
  eta_sec = round(Int, est_total_time - elapsed_time )
  eta = durationstring(eta_sec)
  println("\% Complete: $percentage_complete. ETA: $eta")
  flush(STDOUT)

  Int(any(isnan,sol.u))
end
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptive,1:numRuns)) 
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")



lw = 3

#=
p1 = plot(Δts,fails,ylim=(0,1000),xscale=:log2,yscale=:log10,guidefont=font(16),tickfont=font(16),yguide="Fails Per 1000 Runs",xguide=L"Chosen $\Delta t$",left_margin=100px,top_margin=50px,linewidth=lw,lab=["Euler-Maruyama" "SRIW1"],legendfont=font(14))
p2 = plot(Δts,times,xscale=:log2,yscale=:log10,guidefont=font(16),tickfont=font(14),yguide="Elapsed Time (s)",xguide=L"Chosen $\Delta t$",top_margin=50px,linewidth=lw,lab=["Euler-Maruyama" "SRIW1"],legendfont=font(14))
plot!(Δts,repmat([adaptiveTime],11),linewidth=lw,line=:dash,lab="ESRK+RSwM3")
plot(p1,p2,size=(1200,800))
gui()
=#
p2 = plot(Δts,times,xscale=:log2,yscale=:log2,guidefont=font(16),tickfont=font(14),yguide="Elapsed Time (s)",xguide=L"Chosen $\Delta t$",top_margin=50px,linewidth=lw,lab=["Euler-Maruyama" "SRIW1"],legendfont=font(14))
plot!(Δts,repmat([adaptiveTime],11),linewidth=lw,line=:dash,lab="ESRK+RSwM3",left_margin=75px)
plot(p2,size=(800,800))
gui()
