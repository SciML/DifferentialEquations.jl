######
##FEM Heat Δx Convergence Tests
######
using DifferentialEquations

#Travis CI Test Setting
#Not good plots, but quick for unit tests
Δxs = 1.//2.^(2:-1:1)
Δts = 1//2^(6) * ones(Δxs) #Run at 2^-7 for best plot
#=
# Use this setup to get good plots
Δt = 1//2^(14) #Small Δt for Euler stability, but takes long
N = 4
topΔx = 7
=#

pdeProb = heatProblemExample_moving()

alg = "Euler"; println(alg)
convsim = testConvergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem;alg=alg)

alg = "ImplicitEuler"; println(alg)
convsim2 = testConvergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem;alg=alg)

alg = "CrankNicholson"; println(alg)
convsim3 = testConvergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem;alg=alg)

convplot_fullΔx(convsim,titleStr="")
convplot_fullΔx(convsim2,titleStr="")
convplot_fullΔx(convsim3,titleStr="Dx Convergence Plots")

#Returns true if all converge approximately Δx^2
minimum([convsim.𝒪est["L2"],convsim2.𝒪est["L2"],convsim3.𝒪est["L2"]] - 2 .<.1)
