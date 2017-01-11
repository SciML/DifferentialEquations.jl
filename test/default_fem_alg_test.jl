using DifferentialEquations

### Setup Poisson
dx = 1//2^(3)
mesh = notime_squaremesh([0 1 0 1],dx,:dirichlet)

f = (x) -> sin.(2π.*x[:,1]).*cos.(2π.*x[:,2])
prob = PoissonProblem(f,mesh)

sol = solve(prob)

# Now do Heat

using DiffEqProblemLibrary

#Define a parabolic problem
prob = prob_femheat_birthdeath


#Solve it with a bunch of different algorithms, plot solution
sol = solve(prob)
sol = solve(prob,alg_hints=[:stiff])

true
