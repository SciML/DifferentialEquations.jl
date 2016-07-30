using DifferentialEquations, Plots
srand(100)
prob = linearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,[0,1],Δt=1//2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1,reltol=0)
err1 = sol.errors[:final]
println("Final error for the first solution was $err1")
TEST_PLOT && p1 = plot(sol,plottrue=true)

sol2 =solve(prob::SDEProblem,[0,1],Δt=1//2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1e-1,reltol=0)
err2 = sol2.errors[:final]
println("Final error for the second solution was $err2")
TEST_PLOT && p2 = plot(sol2,plottrue=true)

sol3 =solve(prob::SDEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1e-2,reltol=0)
err3 = sol3.errors[:final]
println("Final error for the third solution was $err3")
TEST_PLOT && p3 = plot(sol3,plottrue=true)

sol4 =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1e-3,reltol=0)
err4 = sol4.errors[:final]
println("Final error for the fourth solution was $err4")
TEST_PLOT && p4 = plot(sol4,plottrue=true)

plot(p1,p2,p3,p4,title="Solutions to Linear SDE at Different Tolerances")
#gui()

true
