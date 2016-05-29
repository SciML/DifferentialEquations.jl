using DifferentialEquations, GrowableArrays
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,1//2^(4),1,fullSave=true,alg="SRI")

fig = PyPlot.figure("pyplot_appx_vs_true",figsize=(10,10))
#First index is the sime, so sol.uFull[1,..] is the initial condition
#Last indices are the indexes of the variables. Since our initial condition
#Has 4 rows and two columns, sol.uFull[..,1] returns the time series for the
#first row, and sol.uFull[..,2] returns the time series for the second.
PyPlot.plot(sol.tFull,sol.uFull[..,1])
PyPlot.plot(sol.tFull,sol.solFull[..,1])
PyPlot.plot(sol.tFull,sol.uFull[..,2])
PyPlot.plot(sol.tFull,sol.solFull[..,2])

## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(14:-1:7) #14->7 good plot

println(@elapsed begin
convsim = testConvergence(Δts,prob,numMonte=Int(5e1),alg="EM")

convsim2 = testConvergence(Δts,prob,numMonte=Int(5e1),alg="RKMil")

convsim3 = testConvergence(Δts,prob,numMonte=Int(5e1),alg="SRI")
end)

abs(convsim.𝒪est["l2"]-.5) + abs(convsim2.𝒪est["l∞"]-1) + abs(convsim3.𝒪est["final"]-1.5)<.2 #High tolerance since low Δts for testing!
