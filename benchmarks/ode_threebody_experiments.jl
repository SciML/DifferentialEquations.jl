using DifferentialEquations
prob = threebodyODEExample()

t₀ = 0.0; T = parse(BigFloat,"17.0652165601579625588917206249")

tspan = [t₀,T]
tspan2 = [t₀,T]
setups = [Dict(:alg=>:DP5)
          Dict(:alg=>:dopri5)
          Dict(:alg=>:ode45)
          Dict(:alg=>:ode45,:norm=>(y)->vecnorm(y,2),:minstep=>1e-40)
          Dict(:alg=>:ExplicitRK)
          Dict(:alg=>:BS5)
          Dict(:alg=>:Tsit5)
          Dict(:alg=>:Feagin14)
          Dict(:alg=>:ode78)
          Dict(:alg=>:DP8)
          Dict(:alg=>:dop853)
          ]

shoot = ode_shootout(prob,tspan,setups,endsol=prob.u₀,abstol=1e-6,reltol=1e-3)

shoot = ode_shootout(prob,tspan,setups,endsol=prob.u₀,abstol=1e-9,reltol=1e-5)

shoot = ode_shootout(prob,tspan2,setups,endsol=prob.u₀,abstol=1e-6,reltol=1e-3)
