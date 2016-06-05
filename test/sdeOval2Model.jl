using DifferentialEquations, Plots, EllipsisNotation, JLD
srand(100)
prob = oval2ModelExample()

sol =solve(prob::SDEProblem,[0;500],Δt=1//2^(10),fullSave=true,alg="SRI",adaptiveAlg="RSwM3",adaptive=true,progressBar=true,saveSteps=1,abstol=1e-6,reltol=1e-4)

#Plots

p1 = plot(sol.tFull,sol.uFull[..,16],top_margin=50px,title="Ecad",xguide="Time",yguide="Concentration",guidefont=font(16),tickfont=font(16))
p2 = plot(sol.tFull,sol.uFull[..,17],top_margin=50px,title="Vim",xguide="Time",yguide="Concentration",guidefont=font(16),tickfont=font(16))
p3 = plot(sol.tFull[1:end-1],diff(sol.tFull)/1000,xguide="Time",yguide="Average Dt over 1000 accepted steps",guidefont=font(16),tickfont=font(16))
plot(p1,p2,p3,layout=@layout([a b;c]),size=(1200,800),title="Adaptive Solution to Stochastic Cell Model")
gui()

save("Oval2Solution.jld","sol",sol,"prob",prob)



##Adaptivity Necessity Tests
sol =solve(prob::SDEProblem,[0;1],Δt=1//2^(8),fullSave=true,alg="EM",adaptive=false,progressBar=true,saveSteps=1,abstol=1e-6,reltol=1e-4)
Int(sol.u[1]!=NaN)
numFails = 0
for i = 1:100
  sol =solve(prob::SDEProblem,[0;1],Δt=1//2^(8),fullSave=true,alg="EM",adaptive=false,progressBar=true,saveSteps=1,abstol=1e-6,reltol=1e-4)
  numFails+=sol.u[1]!=NaN
end
println("The number of Euler-Maruyama Fails is $numFails")
numFails = 0
for i = 1:100
  sol =solve(prob::SDEProblem,[0;1],Δt=1//2^(8),fullSave=true,alg="SRI",adaptive=false,progressBar=true,saveSteps=1,abstol=1e-6,reltol=1e-4)
  numFails+=sol.u[1]!=NaN
end
println("The number of Rossler-SRI Fails is $numFails")
numFails = 0
for i = 1:100
  sol =solve(prob::SDEProblem,[0;1],Δt=1//2^(8),fullSave=true,alg="SRI",adaptiveAlg="RSwM3",adaptive=true,progressBar=true,saveSteps=1,abstol=1e-6,reltol=1e-4)
  numFails+=sol.u[1]!=NaN
end
println("The number of ESRK1+RSwM3 Fails is $numFails")
