@doc """
`HeatProblem`

Wraps the data that define a 2D linear heat equation problem:

``u_t = Δu + f(x,t)``

#Constructors

HeatProblem(sol,Du,f,isLinear): Defines the Dirichlet problem with solution sol,
solution gradient Du = [u_x,u_y], f, and a boolean which states whether the
problem is linear (i.e. linear if f does not depend on u).

HeatProblem(u0,f,gD,gN,isLinear): Defines the problem with initial value u0 (as a function or vector), f,
Dirichlet boundary function gD,  Neumann boundary function gN, and a boolean which states whether the
problem is linear (i.e. linear if f does not depend on u).
""" ->
type HeatProblem <: PdeProblem
  "u0: Initial value function or vector"
  u0
  "Du: Function for the solution gradient [u_x,u_y]"
  Du
  "f: Forcing function in heat equation"
  f
  "gD: Dirichlet boundary data"
  gD
  "gN: Neumann boundary data"
  gN
  "sol: Solution to the heat problem"
  sol
  "knownSol: Boolean which states whether the solution function is given"
  knownSol
  "isLinear: Boolean which states whether the problem is linear or nonlinear"
  isLinear
  σ
  stochastic
  noiseType
  function HeatProblem(sol,Du,f,isLinear;σ=(x)->0,stochastic=false,noiseType="White")
    knownSol = true
    u0(x) = sol(x,0)
    gD = sol
    return(new(u0,Du,f,gD,nothing,sol,knownSol,isLinear,σ,stochastic,noiseType))
  end
  function HeatProblem(u0,f,gD,gN,isLinear;σ=(x)->0,stochastic=false,noiseType="White")
    knownSol = false
    return(new(u0,nothing,f,gD,gN,nothing,knownSol,isLinear,σ,stochastic,noiseType))
  end
end
"""
PoissonProblem

Wraps the data that define a 2D linear Poisson equation problem:

Δu = f(x,t)

#Constructors

PoissonProblem(f,sol,Du,gN,isLinear): Defines the Dirichlet problem with solution sol, solution gradient Du = [u_x,u_y],
f, and Neumann boundary data gN,

PoissonProblem(u0,f,gD,gN,isLinear): Defines the problem with initial value u0 (as a function or vector), f,
Dirichlet boundary function gD, and Neumann boundary function gN.
"""
type PoissonProblem <: PdeProblem
  "f: Forcing function in the Poisson problem"
  f
  "sol: Solution to the Poisson problem"
  sol
  "Du: Gradient of the solution to the Poisson problem"
  Du
  "gD: Dirichlet Boundary Data"
  gD
  "gN: Neumann Boundary Data"
  gN
  "knownSol: Boolean which states whether the solution function is given"
  knownSol
  "isLinear: Boolean which states whether the problem is linear or nonlinear"
  isLinear
  σ
  stochastic
  noiseType
  function PoissonProblem(f,sol,Du,gN,isLinear;σ=(x)->0,stochastic=false,noiseType="White")
    return(new(f,sol,Du,sol,gN,true,isLinear,σ,stochastic,noiseType))
  end
  function PoissonProblem(f,gD,gN,isLinear;σ=(x)->0,stochastic=false,noiseType="White")
    return(new(f,nothing,nothing,gD,gN,false,isLinear,σ,stochastic,noiseType))
  end
end
