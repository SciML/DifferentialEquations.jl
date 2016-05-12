######
##FEM Heat Nonlinear Test
######
using DifferentialEquations

#Define a parabolic problem
T = 2
Δx = 1//2^(4)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
pdeProb = heatProblemExample_birthdeath()

#Solve it with a bunch of different algorithms, plot solution
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
#solplot_appx(res,savefile="plot.svg",zlim=(0,1))
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)
Plots.savefig("plot.pdf")
Plots.savefig("plot.png")
