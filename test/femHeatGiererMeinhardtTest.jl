using DifferentialEquations, Plots

#Define a parabolic problem
T = .25
Δx = 1//2^(4)
Δt = 1//2^(10)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
#prob = heatProblemExample_gierermeinhardt(a=1,α=1,D=[0.01 1.0],ubar=1,vbar=0,β=10) #Saddle
prob = heatProblemExample_gierermeinhardt(D=[9.45e-4 .27],a=10,α=64.5,β=100,ubar=.001) #Spots

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,saveSteps=10)

plot(sol,plottrue=false,zlim=(0,10),cbar=false)
gui()

@linux? animate(sol,zlims=[(0,10),(0,2)],cbar=false) : println("Animation only works with ImageMagick installation, disabled on osx for testing")
