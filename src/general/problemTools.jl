"""
HeatProblem

Wraps the data that define a 2D heat equation problem:

```math
u_t = Δu + f
```

with bounday conditions `gD` on the Dirichlet boundary and gN on the Neumann boundary.
Linearity is determined by whether the forcing function `f` is a function of two
variables (x,t) or three (u,x,t) (with x=[:,1] and y=[:,2]).

If they keyword `σ` is given, then this wraps the data that define a 2D stochastic heat equation

```math
u_t = Δu + f + σdW_t
```

###Constructors

* `HeatProblem(sol,Du,f)`: Defines the Dirichlet problem with solution `sol`,
solution gradient `Du = [u_x,u_y]`, and the forcing function `f`.

* `HeatProblem(u₀,f)`: Defines the problem with initial value `u₀` (as a function) and `f`.
If your initial data is a vector, wrap it as u₀(x) = vector.

Note: If all functions are of (x,t), then the program assumes it's linear. Write
your functions using x = x[:,1] and y = x[:,2].  Use f=f(u,x,t) and σ=σ(u,x,t) (if specified)
for nonlinear problems (with the boundary conditions still (x,t))

###Keyword Arguments

* `gD` = Dirichlet boundary function

* `gN` = Neumann boundary function

* `σ` = The function which multiplies the noise dW. By default σ is 0.

* `noiseType` = A string which specifies the type of noise to be generated. By default
noiseType is "White" for Gaussian Spacetime White Noise.

"""
type HeatProblem <: PdeProblem
  "u₀: Initial value function or vector"
  u₀
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
    u₀(x) = sol(x,0)
    gD = sol
    return(new(u₀,Du,f,gD,gN,sol,knownSol,isLinear,σ,stochastic,noiseType))
  end
  function HeatProblem(u₀,f;gD=(x,t)->zeros(size(x,1)),gN=(x,t)->zeros(size(x,1)),σ=nothing,noiseType="White")
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1))
    else
      stochastic=true
    end
    isLinear = numparameters(f)==2
    knownSol = false
    return(new(u₀,(x)->0,f,gD,gN,(x)->0,knownSol,isLinear,σ,stochastic,noiseType))
  end
end

doc"""
PoissonProblem

Wraps the data that define a 2D linear Poisson equation problem:

```math
-Δu = f
```

with bounday conditions `gD` on the Dirichlet boundary and gN on the Neumann boundary.
Linearity is determined by whether the forcing function `f` is a function of two
variables (x,t) or three (u,x,t) (with x=[:,1] and y=[:,2]).

If they keyword `σ` is given, then this wraps the data that define a 2D stochastic heat equation

```math
-Δu = f + σdW
```

###Constructors

PoissonProblem(f,sol,Du): Defines the Dirichlet problem with solution `sol`, solution gradient `Du = [u_x,u_y]`,
and forcing function `f`

PoissonProblem(u₀,f): Defines the problem with initial value `u₀` (as a function) and f.
If your initial data is a vector, wrap it as u₀(x) = vector.

Note: If all functions are of (x,t), then the program assumes it's linear. Write
your functions using x = x[:,1] and y = x[:,2].  Use f=f(u,x,t) and σ=σ(u,x,t) (if specified)
for nonlinear problems (with the boundary conditions still (x,t))

###Keyword Arguments

* `gD` = Dirichlet boundary function

* `gN` = Neumann boundary function

* `σ` = The function which multiplies the noise ``dW``. By default `σ` is 0.

* `noiseType` = A string which specifies the type of noise to be generated. By default
`noiseType` is "White" for Gaussian Spacetime White Noise.

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

"""
numparameters(f)

Returns the number of parameters of `f` for the method which has the most parameters.
"""
function numparameters(f)
  if length(methods(f))>1
    warn("Number of methods for f is greater than 1. Choosing linearity based off of method with most parameters")
  end
  maximum([length(m.sig.parameters) for m in methods(f)])
end
