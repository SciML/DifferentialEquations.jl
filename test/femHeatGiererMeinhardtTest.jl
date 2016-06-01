using DifferentialEquations, Plots

#Define a parabolic problem
T = 100
Δx = 1//2^(4)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_gierermeinhardt()

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,saveSteps=100)

plot(sol,plottrue=false,zlim=(0,20),cbar=false)
gui()

animate(sol,zlim=(0,20),cbar=false)
