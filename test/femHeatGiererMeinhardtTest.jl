using DifferentialEquations, Plots

#Define a parabolic problem
T = .5
k = 4
Δx = 1//2^(k)
Δt = 1//2^(10)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
#prob = heatProblemExample_gierermeinhardt(a=1,α=1,D=[0.01 1.0],ubar=1,vbar=0,β=10) #Saddle
prob = heatProblemExample_gierermeinhardt(D=[9.45e-4 .27],a=10,α=64.5,β=100,ubar=.001) #Spots

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,saveSteps=10)

plot(sol,plottrue=false,zlim=(0,10),cbar=false)
gui()

animate(sol,zlims=[(0,10),(0,2)],cbar=false)

Plots.heatmap(vec(reshape(sol.femMesh.node[:,1],2^(4)+1,2^(k)+1)[1,:]),
              vec(reshape(sol.femMesh.node[:,2],2^(4)+1,2^(k)+1)[:,1]),
              reshape(sol.timeSeries[1][end],2^(4)+1,2^(k)+1))
zlims=[(0,10),(0,2)]
cbar = false

atomLoaded = isdefined(Main,:Atom)
Plots.pyplot(reuse=true,size=(1500,750))
@gif for j=1:length(sol.timeSeries[1])
  ps = Any[]
  for i=1:sol.prob.numVars
    push!(ps,Plots.heatmap(vec(reshape(sol.femMesh.node[:,1],2^(4)+1,2^(k)+1)[1,:]),
                  vec(reshape(sol.femMesh.node[:,2],2^(4)+1,2^(k)+1)[:,1]),
                  reshape(sol.timeSeries[i][j],2^(4)+1,2^(k)+1),zlim=zlims[i],cbar=cbar))
  end
  plot(ps...)
  atomLoaded ? Main.Atom.progress(j/length(sol.timeSeries[1])) : nothing #Use Atom's progressbar if loaded
end
