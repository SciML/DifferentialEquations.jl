using DifferentialEquations, SIUnits, SIUnits.ShortUnits
f = (t,y) -> 0.5*y*t / 3.0s
u₀ = 1.5Newton
prob = ODEProblem(f,u₀)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:Midpoint)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:ExplicitRK,adaptive=true)

plot(sol)
#plot(sol.t,sol[:]) # Doesn't work because no unit recipe
#plot(sol) #Doesn't work because the one above doesn't work

β = 0.6
σ(y,t) = β*y
prob = SDEProblem(f,σ,u₀)

sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:EM)

### Setup
Δx = (1//2^(3))Meter
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet) #Fails at meshgrid because unit ranges not supported

#=
Δx = (1//2^(3))
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet) #Fails at meshgrid because unit ranges not supported

f(x)=sin(2π.*x[:,1]).*cos(2π.*x[:,2])Meter
prob = PoissonProblem(f)

sol = solve(fem_mesh,prob)
=#

#Define a parabolic problem
T = 1
Δx = 1//2^(3)
Δt = 1//2^(7)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:dirichlet)
f(u,x,t)  = (float(ones(size(x,1))))Newton - .5u
u₀(x) = map((x)->(x)Newton,zeros(size(x,1)))
prob = HeatProblem(u₀,f)
#Solve it with a bunch of different algorithms, plot solution
println("Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler)
