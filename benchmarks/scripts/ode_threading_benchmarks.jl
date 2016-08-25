using DifferentialEquations
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = 1.01*u[i]
  end
end
analytic = (t,u₀) -> u₀*exp(1.01*t)
problarge = ODEProblem(f,rand(100,100),analytic=analytic)



elapsed1 = @elapsed sol1 =solve(problarge::ODEProblem,[0,10];alg=:DP5)

elapsed1 = @elapsed sol1 =solve(problarge::ODEProblem,[0,10];alg=:DP5Threaded)


setups = [Dict(:alg=>:DP5)
          Dict(:alg=>:DP5Threaded)
          Dict(:abstol=>1e-3,:reltol=>1e-6,:alg=>:ode45) # Fix ODE to be normal
          Dict(:alg=>:dopri5)]
]
names = ["DifferentialEquations";"DEThreaded";"ODE";"ODEInterface"]
shoot = ode_shootout(prob,tspan,setups;Δt=1/2^(10),names=names)
