using DifferentialEquations, Plots
srand(100)
prob = prob_sde_2Dlinear

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1,reltol=0)
err1 = sol.errors[:final]
TEST_PLOT && p1 = plot(sol,plot_analytic=true,legend=false,title="tol = 1")

println("1e-1")
sol2 =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1e-1,reltol=0)
err2 = sol2.errors[:final]
TEST_PLOT && p2 = plot(sol2,plot_analytic=true,legend=false,title="tol = 1e-1")

println("1e-2")
sol3 =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1e-2,reltol=0)
err3 = sol3.errors[:final]
TEST_PLOT && p3 = plot(sol3,plot_analytic=true,legend=false,title="tol = 1e-2")

println("1e-3")
sol4 =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1e-3,reltol=0)
err4 = sol4.errors[:final]
TEST_PLOT && p4 = plot(sol4,plot_analytic=true,legend=false,title="tol = 1e-3")

TEST_PLOT && plot(p1,p2,p3,p4,title="Solutions to Linear SDE at Different Tolerances",size=(1200,800))
#gui()

println("""
Final error for the solutions were:
          $err1
          $err2
          $err3
          $err4""")

err4 < err2 && err3 < err1
