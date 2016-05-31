using DifferentialEquations, Plots

#Define a parabolic problem
T = 40000
Δx = 1//2^(3)
Δt = 1//2^(9)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_gierermeinhardt()

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,saveSteps=1000)

plot(sol,plottrue=false,zlim=(0,20),cbar=false)
gui()
