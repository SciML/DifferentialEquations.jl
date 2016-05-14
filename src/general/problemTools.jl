@doc_str"""
HeatProblem

Wraps the data that define a 2D linear heat equation problem:

$$u_t = Δu + f(x,t)$$

#Constructors

HeatProblem(sol,Du,f,isLinear): Defines the Dirichlet problem with solution sol,
solution gradient Du = [u_x,u_y], f, and a boolean which states whether the
problem is linear (i.e. linear if f does not depend on u).

HeatProblem(u0,f,gD,gN,isLinear): Defines the problem with initial value u0 (as a function or vector), f,
Dirichlet boundary function gD,  Neumann boundary function gN, and a boolean which states whether the
problem is linear (i.e. linear if f does not depend on u).

Note: If isLinear is true, then all functions must only be functions of (x,t). If
isLinear is false, then f=f(u,x,t) and σ=σ(u,x,t) (if specified), while the other
functions are only functions of (x,t).

#Keyword Arguments

The constructors take the following keyword arguments:

σ = The function which multiplies the noise dW. By default σ is 0.

stochastic = A boolean which specifies if the problem is stochastic. By default
stochastic is false.

noiseType = A string which specifies the type of noise to be generated. By default
noiseType is "White" for Gaussian Spacetime White Noise.

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
  function HeatProblem(sol,Du,f;gN=(x,t)->zeros(size(x,1)),σ=nothing,noiseType="White")
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1))
    else
      stochastic=true
    end
    isLinear = numparameters(f)==2
    knownSol = true
    u0(x) = sol(x,0)
    gD = sol
    return(new(u0,Du,f,gD,gN,sol,knownSol,isLinear,σ,stochastic,noiseType))
  end
  function HeatProblem(u0,f;gD=(x,t)->zeros(size(x,1)),gN=(x,t)->zeros(size(x,1)),σ=nothing,noiseType="White")
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1))
    else
      stochastic=true
    end
    isLinear = numparameters(f)==2
    knownSol = false
    return(new(u0,(x)->0,f,gD,gN,(x)->0,knownSol,isLinear,σ,stochastic,noiseType))
  end
end
@doc_str"""
PoissonProblem

Wraps the data that define a 2D linear Poisson equation problem:

``$$#Δu = f(x,t)$$``

#Constructors

PoissonProblem(f,sol,Du,gN,isLinear): Defines the Dirichlet problem with solution sol, solution gradient Du = [u_x,u_y],
f, and Neumann boundary data gN,

PoissonProblem(u0,f,gD,gN,isLinear): Defines the problem with initial value u0 (as a function or vector), f,
Dirichlet boundary function gD, and Neumann boundary function gN.

Note: If isLinear is true, then all functions must only be functions of (x). If
isLinear is false, then f=f(u,x) and σ=σ(u,x) (if specified), while the other
functions are only functions of (x).

#Keyword Arguments

`σ` = The function which multiplies the noise ``dW``. By default `σ` is 0.

noiseType = A string which specifies the type of noise to be generated. By default
noiseType is "White" for Gaussian Spacetime White Noise.

"""
type PoissonProblem <: PdeProblem
  "f: Forcing function in the Poisson problem"
  f::Function
  "sol: Solution to the Poisson problem"
  sol::Function
  "Du: Gradient of the solution to the Poisson problem"
  Du::Function
  "gD: Dirichlet Boundary Data"
  gD::Function
  "gN: Neumann Boundary Data"
  gN::Function
  "knownSol: Boolean which states whether the solution function is given"
  knownSol::Bool
  "isLinear: Boolean which states whether the problem is linear or nonlinear"
  isLinear::Bool
  σ::Function
  stochastic::Bool
  noiseType::AbstractString
  function PoissonProblem(f,sol,Du;gN=(x)->zeros(size(x,1)),σ=nothing,noiseType="White")
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1))
    else
      stochastic=true
    end
    isLinear = numparameters(f)==1
    return(new(f,sol,Du,sol,gN,true,isLinear,σ,stochastic,noiseType))
  end
  function PoissonProblem(f;gD=(x)->zeros(size(x,1)),gN=(x)->zeros(size(x,1)),σ=nothing,noiseType="White")
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1))
    else
      stochastic = true
    end

    isLinear = numparameters(f)==1
    return(new(f,(x)->0,(x)->0,gD,gN,false,isLinear,σ,stochastic,noiseType))
  end
end

function numparameters(f)
  if length(methods(f))>1
    warn("Number of methods for f is greater than 1. Choosing linearity based off of method with most parameters")
  end
  maximum([length(m.sig.parameters) for m in methods(f)])
end
