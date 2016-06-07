using DifferentialEquations, Plots
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=1,reltol=0)
err1 = sol.errors["final"]
p1 = plot(sol,plottrue=true,legend=false,title="tol = 1")

println("1e-1")
sol2 =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=1e-1,reltol=0)
err2 = sol2.errors["final"]
p2 = plot(sol2,plottrue=true,legend=false,title="tol = 1e-1")

println("1e-2")
sol3 =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=1e-2,reltol=0)
err3 = sol3.errors["final"]
p3 = plot(sol3,plottrue=true,legend=false,title="tol = 1e-2")

println("1e-3")
sol4 =solve(prob::SDEProblem,[0,1],Δt=1/2^(4),fullSave=true,alg="SRI",adaptive=true,abstol=1e-3,reltol=0)
err4 = sol4.errors["final"]
p4 = plot(sol4,plottrue=true,legend=false,title="tol = 1e-3")

plot(p1,p2,p3,p4,title="Solutions to Linear SDE at Different Tolerances",size=(1200,800))
#gui()

println("""
Final error for the solutions were:
          $err1
          $err2
          $err3
          $err4""")

err4 < err2 && err3 < err1
