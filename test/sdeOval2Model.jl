using DifferentialEquations, Plots, EllipsisNotation, LaTeXStrings
srand(100)
prob = oval2ModelExample(largeFluctuations=true,useBigs=false)


## Big Run
@time sol = solve(prob::SDEProblem,[0;500],Δt=(1/2)^(8),fullSave=true,alg="SRIW1Optimized",
          adaptiveAlg="RSwM3",adaptive=true,progressBar=true,
          saveSteps=1000,abstol=1e-5,reltol=1e-3)

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

#using JLD
#save("Oval2Solution.jld","sol",sol,"prob",prob)
