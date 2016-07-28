######
##FEM Heat Δx Convergence Tests
######
using DifferentialEquations, Plots

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

prob = heatProblemExample_moving()

alg=:Euler; println(alg)
sim = test_convergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem,Δxs;alg=alg)

alg=:ImplicitEuler; println(alg)
sim2 = test_convergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem,Δxs;alg=alg)

alg=:CrankNicholson; println(alg)
sim3 = test_convergence(Δts::AbstractArray,Δxs::AbstractArray,prob::HeatProblem,Δxs;alg=alg)

plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,400))

#Returns true if all converge approximately Δx^2
minimum([sim.𝒪est[:L2],sim2.𝒪est[:L2],sim3.𝒪est[:L2]] - 2 .<.1)
