######
##FEM Heat Δt Convergence Tests
######
using DifferentialEquations, Plots

#Convergences estimate has not converged in this range
#Should decrease Δx/Δt for better estimate
N = 2 #Number of different Δt to solve at, 2 for test speed
topΔt = 6 # 1//2^(topΔt-1) is the max Δt. Small for test speed
prob = heatProblemExample_moving() #also try heatProblemExample_pure() or heatProblemExample_diffuse()
Δts = 1.//2.^(topΔt-1:-1:N)
Δxs = 1//2^(5) * ones(Δts) #Run at 2^-7 for best plot


alg = "Euler"; println(alg) #Unstable due to μ
sim = testConvergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem,Δts;alg=alg)

alg = "ImplicitEuler"; println(alg)
sim2 = testConvergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem,Δts;alg=alg)

alg = "CrankNicholson"; println(alg) #Bound by spatial discretization error at low Δt, decrease Δx for full convergence
sim3 = testConvergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem,Δts;alg=alg)

plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,400))
#Note: Stabilizes in H1 due to high Δx-error, reduce Δx and it converges further.

#Returns true if ImplicitEuler converges like Δt and
#CN convergeces like >Δt^2 (approaches Δt^2 as Δt and Δx is smaller
minimum([abs(sim2.𝒪est["L2"]-1)<.3 sim3.𝒪est["L2"]>2])
