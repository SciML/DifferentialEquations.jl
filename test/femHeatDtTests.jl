######
##FEM Heat Δt Convergence Tests
######
using DifferentialEquations

T = 1
Δx = 1//2^(3) #Run at 2^-6 for better plot, at 2^-5 for test times
N = 4
topΔt = 8
pdeProb = heatProblemExample_moving() #also try heatProblemExample_pure() or heatProblemExample_diffuse()


alg = "Euler" #Unstable due to μ
solutions = cell(N)
for i = 1:N
  Δt = 1//2^(topΔt-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg)
  solutions[i] = res
end
simres = ConvergenceSimulation(solutions)
if !isdefined(:testState) #Don't plot during test
  convplot_fullΔt(simres,titleStr="Euler Convergence Plots",savefile="eulerdtconv.svg")
end

alg = "ImplicitEuler"
solutions = cell(N)
for i = 1:N
  Δt = 1//2^(topΔt-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg,solver="LU")
  solutions[i] = res
end
simres2 = ConvergenceSimulation(solutions)
if !isdefined(:testState) #Don't plot during test
  convplot_fullΔt(simres2,titleStr="Implicit Euler Convergence Plots",savefile="impeulerdtconv.svg")
end

alg = "CrankNicholson" #Bound by spatial discretization error at low Δt, decrease Δx for full convergence
solutions = cell(N)
for i = 1:N
  Δt = 1//2^(topΔt-i)
  femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
  res = fem_solveheat(femMesh::FEMmesh,pdeProb,alg=alg,solver="GMRES") #LU faster, but GMRES more stable
  solutions[i] = res
end
simres3 = ConvergenceSimulation(solutions)
if !isdefined(:testState) #Don't plot during test
  convplot_fullΔt(simres3,titleStr="Crank-Nicholson Convergence Plots",savefile="crankdtconv.svg")
end
#Note: Stabilizes in H1 due to high Δx-error, reduce Δx and it converges further.
