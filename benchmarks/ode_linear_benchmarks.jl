using DifferentialEquations
probnum = linearODEExample()
prob = twoDimlinearODEExample!(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
tspan = [0,1]
using BenchmarkTools

setups = [Dict(:alg=>:RK4);Dict(:alg=>:ode4)]

shoot = ode_shootout(probnum,tspan,setups;Δt=1/2^(6))
shoot = ode_shootout(prob,tspan,setups;Δt=1/2^(6))

## Standard Tolerance
tspan = [0,10]
setups = [Dict(:alg=>:DP5)
          Dict(:abstol=>1e-3,:reltol=>1e-6,:alg=>:ode45) # Fix ODE to be normal
          Dict(:alg=>:dopri5)]
names = ["DifferentialEquations";"ODE";"ODEInterface"]
shoot = ode_shootout(prob,tspan,setups;Δt=1/2^(10),names=names)

# Low Tolerance

setups = [Dict(:alg=>:DP5)
          Dict(:alg=>:ode45)
          Dict(:alg=>:dopri5)]

shoot = ode_shootout(prob,tspan,setups;Δt=1/2^(10),names=names,abstol=1e-6,roltol=1e-6)

# Number

shoot = ode_shootout(probnum,tspan,setups;Δt=1/2^(10),names=names,abstol=1e-6,roltol=1e-6)

## Test how much choices matter

setups = [Dict(:alg=>:DP5)
          Dict(:alg=>:DP5Vectorized)
          Dict(:alg=>:DP5,:timechoicealg=>:Simple)
          Dict(:alg=>:DP5,:fullnormalize=>true)
          Dict(:alg=>:ExplicitRK)
          Dict(:alg=>:DP5,:β=>0.04)
          Dict(:alg=>:BS5)]

names = ["Standard DP5","Vectorized","No PI Control","Full Normalized","Tableau Form","Lower β","BS5"]

shoot = ode_shootout(prob,tspan,setups,names=names)

## Test at ODE.jl defaults

setups = [Dict(:abstol=>1e-8,:reltol=>1e-5,:alg=>:DP5)
          Dict(:alg=>:ode45)]
names = ["DifferentialEquations","ODE"]
shoot = ode_shootout(prob,tspan,setups;names=names)

# Other Methods

setups = [Dict(:alg=>:DP5)
          Dict(:alg=>:BS3)
          Dict(:alg=>:BS5)
          Dict(:alg=>:Tsit5)
          Dict(:alg=>:DP8)
          Dict(:alg=>:dop853)
          Dict(:alg=>:ode78)]

shoot = ode_shootout(prob,tspan,setups)
