using DifferentialEquations, GrowableArrays
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,1//2^(4),1,fullSave=true,alg="SRI")
monteCarloSim(1//2^(4),prob::SDEProblem,T=1)
#First index is the sime, so sol.uFull[1,..] is the initial condition
#Last indices are the indexes of the variables. Since our initial condition
#Has 4 rows and two columns, sol.uFull[..,1] returns the time series for the
#first row, and sol.uFull[..,2] returns the time series for the second.
plot(sol,plottrue=true)
# reshape(sol.uFull,size(sol.uFull,1),size(sol.uFull,2)*size(sol.uFull,3))
PyPlot.plot(sol.tFull,sol.uFull[..,1])
PyPlot.plot(sol.tFull,sol.solFull[..,1])
PyPlot.plot(sol.tFull,sol.uFull[..,2])
PyPlot.plot(sol.tFull,sol.solFull[..,2])

## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(14:-1:7) #14->7 good plot

println(@elapsed begin
sim = testConvergence(Δts,prob,numMonte=Int(5e1),alg="EM")

sim2 = testConvergence(Δts,prob,numMonte=Int(5e1),alg="RKMil")

sim3 = testConvergence(Δts,prob,numMonte=Int(5e1),alg="SRI")
end)

abs(sim.𝒪est["l2"]-.5) + abs(sim2.𝒪est["l∞"]-1) + abs(sim3.𝒪est["final"]-1.5)<.2 #High tolerance since low Δts for testing!
