######
##FEM Heat Δt Convergence Tests
######
using DifferentialEquations

#Convergences estimate has not converged in this range
#Should decrease Δx/Δt for better estimate
T = 1
Δx = 1//2^(5) #Run at 2^-7 for best plot
N = 2 #Number of different Δt to solve at, 2 for test speed
topΔt = 5 # 1//2^(topΔt-1) is the max Δt. Small for test speed
pdeProb = heatProblemExample_moving() #also try heatProblemExample_pure() or heatProblemExample_diffuse()

alg = "Euler" #Unstable due to μ
println(alg)
solutions = cell(N)
for i = 1:N
  Δt = 1//2^(topΔt-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg)
  solutions[i] = res
end
simres = ConvergenceSimulation(solutions)

alg = "ImplicitEuler"
println(alg)
solutions = cell(N)
for i = 1:N
  Δt = 1//2^(topΔt-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg,solver="LU")
  solutions[i] = res
end
simres2 = ConvergenceSimulation(solutions)

alg = "CrankNicholson" #Bound by spatial discretization error at low Δt, decrease Δx for full convergence
println(alg)
solutions = cell(N)
for i = 1:N
  Δt = 1//2^(topΔt-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg,solver="LU") #LU faster, but GMRES more stable
  solutions[i] = res
end
simres3 = ConvergenceSimulation(solutions)

convplot_fullΔt(simres3,titleStr="Crank-Nicholson Convergence Plots")
#Note: Stabilizes in H1 due to high Δx-error, reduce Δx and it converges further.

#Returns true if ImplicitEuler converges like Δt and
#CN convergeces like >Δt^2 (approaches Δt^2 as Δt and Δx is smaller
minimum([abs(simres2.ConvEst_l2-1)<.3 simres3.ConvEst_l2>2])
