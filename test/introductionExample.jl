######
##FEM Poisson Method Tests
######
using DifferentialEquations

Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")
pdeProb = poissonProblemExample_wave()

sol = solve(femMesh,pdeProb,solver="CG")

plot(sol::FEMSolution,plottrue=false) #To save the plot, use savefig("plot.png") or "plot.pdf", etc.

#Error should be low
sol.errors["L2"] < 1e-4
