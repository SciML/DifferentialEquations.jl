using DifferentialEquations, Plots, EllipsisNotation, LaTeXStrings
srand(100)
prob = oval2ModelExample(largeFluctuations=true,useBigs=false)

## Big Run
@time sol = solve(prob::SDEProblem,[0;500],Δt=(1/2)^(8),fullSave=true,alg=:SRIW1Optimized,
          adaptivealg=:RSwM3,adaptive=true,progressBar=true,qmax=4,
          saveSteps=1000,abstol=1e-5,reltol=1e-3)

#Plots
lw = 2
lw2 = 3

js = 16:21
Δts = 1./2.^(js)

adaptiveTime = 429.167375442

fails = [132 131 78
          30 26 17
          5  6  0
          1  1  0
          0  0  0
          0  0  0]

times = [133.346231197  211.915111975  609.267414537
         269.089020237  428.278565948  1244.056558928
         580.135419189  861.005910096  2491.373020819
         1138.412863403 1727.91452979  4932.702912547
         2286.353072179 3439.895966342 9827.155740715
         4562.201874801 6891.354640001 19564.165421295]


p1 = plot(sol.tFull,sol.uFull[..,16],top_margin=50px,title="(A) Timeseries of Ecad Concentration",xguide="Time (s)",yguide="Concentration",guidefont=font(16),tickfont=font(16),linewidth=lw,left_margin=85px,leg=false)
p2 = plot(sol.tFull,sol.uFull[..,17],top_margin=50px,title="(B) Timeseries of Vim Concentration",xguide="Time (s)",yguide="Concentration",guidefont=font(16),tickfont=font(16),linewidth=lw,leg=false)
p3 = plot(sol.tFull,sol.ΔtFull,xguide="Time (s)",yguide="Accepted Dt",title="(C) Accepted Dt vs Time", guidefont=font(16),tickfont=font(16),yscale=:log10,linewidth=lw,left_margin=110px,bottom_margin=65px,leg=false)
p4 =  plot(Δts,times,xscale=:log2,yscale=:log2,guidefont=font(16),tickfont=font(14),yguide="Elapsed Time (s)",xguide=L"Chosen $\Delta t$",linewidth=lw2,lab=["Euler-Maruyama" "RK-Mil" "RosslerSRI"],legendfont=font(12),title="(D) Fixed Timestep vs Adaptive Timing Test")
plot!(Δts,repmat([adaptiveTime],11),linewidth=lw2,line=:dash,lab="ESRK+RSwM3",left_margin=75px)
scatter!([2.0^(-20);2.0^(-20);2.0^(-18)],[times[5,1];times[5,2];times[3,3]],markersize=20,c=:red,lab="")

plot(p1,p2,p3,p4,layout=@layout([a b;c d]),size=(1200,800))
gui()
