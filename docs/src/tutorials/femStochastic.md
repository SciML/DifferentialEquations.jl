# Stochastic Finite Element Examples

### Finite Element Stochastic Poisson Equation

We can solve the same PDE as a stochastic PDE, -Δu=f+gdW, with additive space-time white noise by specifying the problem as:

```julia
"Example problem with deterministic solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_noisyWave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  gN(x) = 0
  isLinear = true
  stochastic = true
  σ(x) = 5 #Additive noise, a big amount!
  return(PoissonProblem(f,sol,Du,gN,isLinear,σ=σ,stochastic=stochastic))
end
```

and using the same solving commands as shown in [femStochasticPoissonSolo.jl](/src/femStochasticPoissonSolo.jl). This gives the following plot:

<img src="/src/examples/introductionStochasticExample.png" width="750" align="middle" />

### Finite Element Stochastic Heat Equation

The last equation we will solve in this introductory example will be a nonlinear stochastic heat equation u_t=Δu+f+gdW with forcing function `f(u)=.5-u`, noise function `g(u)=100u^2` and
initial condition `u0=0`. We would expect this system to rise towards the deterministic steady state `u=2` (but stay in mean a bit below it due to 1st order "Milstein" effects), gaining more noise as it increases. This is specified as follows:

```julia
"Example problem which starts with 0 and solves with f(u)=1-.1u"
function heatProblemExample_stochasticbirthdeath()
  gD(x,t) = zeros(size(x,1))
  f(u,x,t)  = ones(size(x,1)) - .5u
  u0(x) = zeros(size(x,1))
  gN(x,t) = 0
  isLinear = false
  stochastic = true
  σ(u,x,t) = 100u.^2
  return(HeatProblem(u0,f,gD,gN,isLinear,σ=σ,stochastic=stochastic))
end
```

As shown in [femStochasticHeatAnimationTest.jl](/src/femStochasticHeatAnimationTest.jl), we use the following code create an animation of the solution:

```julia
T = 5
Δx = 1//2^(4)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
pdeProb = heatProblemExample_stochasticbirthdeath()

res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler",fullSave=true)
solplot_animation(res::FEMSolution;zlim=(0,2),vmax=.1,cbar=false)
```

<img src="/src/examples/stochasticHeatAnimation.gif" width="750" align="middle" />
