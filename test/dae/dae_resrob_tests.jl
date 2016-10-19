using DifferentialEquations

prob = prob_dae_resrob
tspan = [0;1]
sol = solve(prob,tspan)
#plot(sol)

true
