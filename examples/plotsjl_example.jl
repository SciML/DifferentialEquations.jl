using DifferentialEquations, Plots
srand(10)

prob = vanDerPolExample()
sol = solve(prob,[0,20])
p1 = plot(sol,title="Solution to the Van der Pol Equations",lw=3,leg=false)

prob = twoDimlinearSDEExample()

sol = solve(prob::SDEProblem,[0,1],Δt=1/2^(8),save_timeseries=true,alg=:SRI)
p2 = plot(sol,leg=false,lw=2,title="8 Independent Black-Scholes Solutions",xaxis ="Time",yaxis="Price")

Δx = 1//2^(4)
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet)

f(x)=sin(2π.*x[:,1]).*cos(2π.*x[:,2])
prob = PoissonProblem(f)

sol = solve(fem_mesh,prob)

p3 = plot(sol,title="Finite-Element Solution to the Poisson Equation",cbar=false)

plot(p1,p3,p2,layout=@layout([a b;c]),size=(1000,800))
gui()
