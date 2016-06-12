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

* `numVars` = Number of variables in the system. Automatically calculated from u₀ in most cases.

* `D` = Array which defines the diffusion coefficients. Defaults to 1's.
"""
type HeatProblem <: DEProblem
  "u₀: Initial value function"
  u₀::Function
  "Du: Function for the solution gradient [u_x,u_y]"
  Du::Function
  "f: Forcing function in heat equation"
  f::Function
  "gD: Dirichlet boundary data"
  gD#::Function
  "gN: Neumann boundary data"
  gN#::Function
  "sol: Solution to the heat problem"
  sol::Function
  "knownSol: Boolean which states whether the solution function is given"
  knownSol::Bool
  "isLinear: Boolean which states whether the problem is linear or nonlinear"
  isLinear::Bool
  numVars::Int
  σ::Function
  stochastic::Bool
  noiseType::AbstractString
  D#AbstractArray
  function HeatProblem(sol,Du,f;gN=nothing,σ=nothing,noiseType="White",numVars=nothing,D=nothing)
    isLinear = numparameters(f)==2
    knownSol = true
    u₀(x) = sol(x,0)
    numVars = size(u₀([0 0
                       0 0
                       0 0]),2)
    gD = sol
    if gN == nothing
      gN=(x,t)->zeros(size(x,1),numVars)
    end
    if σ==nothing
      stochastic=false
      σ=(x,t)->zeros(size(x,1),numVars)
    else
      stochastic=true
    end
    if D == nothing
      if numVars == 1
        D = 1.0
      else
        D = ones(1,numVars)
      end
    end
    return(new(u₀,Du,f,gD,gN,sol,knownSol,isLinear,numVars,σ,stochastic,noiseType,D))
  end
  function HeatProblem(u₀,f;gD=nothing,gN=nothing,σ=nothing,noiseType="White",numVars=nothing,D=nothing)
    if σ==nothing
      stochastic=false
      σ=(x,t)->zeros(size(x,1))
    else
      stochastic=true
    end
    isLinear = numparameters(f)==2
    knownSol = false
    if isLinear
      if u₀==nothing
        u₀=(x)->zeros(size(x,1))
      end
      if gD == nothing
        gD=(x,t)->zeros(size(x,1))
      end
      if gN == nothing
        gN=(x,t)->zeros(size(x,1))
      end
      if D == nothing
        D = 1.0
      end
      numVars = 1
    end
    if !isLinear #nonlinear
      if u₀==nothing && numVars == nothing
        warn("u₀ and numVars must be given. numVars assumed 1.")
        numVars = 1
        u₀=(x)->zeros(size(x,1),numVars)
        if gD == nothing
          gD=(x,t)->zeros(size(x,1),numVars)
        end
        if gN == nothing
          gN=(x,t)->zeros(size(x,1),numVars)
        end
        if D == nothing
          D = 1.0
        end
      elseif u₀==nothing #numVars!=nothing
        u₀=(x)->zeros(size(x,1),numVars) #Default to zero
        if gD == nothing
          gD=(x,t)->zeros(size(x,1),numVars)
        end
        if gN == nothing
          gN=(x,t)->zeros(size(x,1),numVars)
        end
        if D == nothing
          D = ones(1,numVars)
        end
      elseif numVars==nothing #If u₀ is given but numVars is not, we're still okay. Generate from size in function.
        numVars=0 #Placeholder, update gD and gN in solver
      end
    end
    return(new(u₀,(x)->0,f,gD,gN,(x)->0,knownSol,isLinear,numVars,σ,stochastic,noiseType,D))
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

* `numVars` = The number of variables in the Poisson system. Automatically calculated in many cases.

* `D` = Vector of diffusion coefficients. Defaults to ones.

"""
type PoissonProblem <: DEProblem
  "f: Forcing function in the Poisson problem"
  f::Function
  "sol: Solution to the Poisson problem"
  sol::Function
  "Du: Gradient of the solution to the Poisson problem"
  Du::Function
  "gD: Dirichlet Boundary Data"
  gD#::Nullable{Function}
  "gN: Neumann Boundary Data"
  gN#::Nullable{Function}
  "knownSol: Boolean which states whether the solution function is given"
  knownSol::Bool
  "isLinear: Boolean which states whether the problem is linear or nonlinear"
  isLinear::Bool
  u₀::Function
  numVars::Int
  σ::Function
  stochastic::Bool
  noiseType::AbstractString
  D#::AbstractArray
  function PoissonProblem(f,sol,Du;gN=nothing,σ=nothing,u₀=nothing,noiseType="White",numVars=nothing,D=nothing)
    gD = sol
    numVars = size(sol([0 0
                        0 0
                        0 0]),2)
    isLinear = numparameters(f)==1
    if gN == nothing
      gN=(x)->zeros(size(x,1),numVars)
    end
    if u₀==nothing
      u₀=(x)->zeros(size(x,1),numVars)
    end
    if D == nothing
      if numVars == 1
        D = 1.0
      else
        D = ones(1,numVars)
      end
    end
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1),numVars)
    else
      stochastic=true
    end
    return(new(f,sol,Du,sol,gN,true,isLinear,u₀,numVars,σ,stochastic,noiseType,D))
  end
  function PoissonProblem(f;gD=nothing,gN=nothing,u₀=nothing,σ=nothing,noiseType="White",numVars=nothing,D=nothing)
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1))
    else
      stochastic = true
    end
    isLinear = numparameters(f)==1
    if isLinear && u₀==nothing
      u₀=(x)->zeros(size(x,1))
      if gD == nothing
        gD=(x)->zeros(size(x,1))
      end
      if gN == nothing
        gN=(x)->zeros(size(x,1))
      end
      if D == nothing
        D = 1.0
      end
      numVars = 1
    end
    if !isLinear #nonlinear
      if u₀==nothing && numVars == nothing
        warn("u₀ and numVars must be given. numVars assumed 1.")
        numVars = 1
        u₀=(x)->zeros(size(x,1))
        if gD == nothing
          gD=(x)->zeros(size(x,1))
        end
        if gN == nothing
          gN=(x)->zeros(size(x,1))
        end
        if D == nothing
          D = 1.0
        end
      elseif u₀==nothing #numVars!=nothing
        u₀=(x)->zeros(size(x,1),numVars) #Default to zero
        if gD == nothing
          gD=(x)->zeros(size(x,1),numVars)
        end
        if gN == nothing
          gN=(x)->zeros(size(x,1),numVars)
        end
        if D == nothing
          D = ones(1,numVars)
        end
      elseif numVars==nothing #If u₀ is given but numVars is not, we're still okay. Generate from size in function.
        numVars=0 #Placeholder, update gD and gN in solver
      end
    end
    return(new(f,(x)->0,(x)->0,gD,gN,false,isLinear,u₀,numVars,σ,stochastic,noiseType,D))
  end
end

"""
SDEProblem

Wraps the data which defines an SDE problem

```math
u = f(u,t)dt + Σσᵢ(u,t)dWⁱ
```

with initial condition u₀.

### Fields

* `f`: The drift function in the SDE.
* `σ`: The noise function in the SDE.
* `u₀`: The initial condition.
* `sol`: A function which describes the solution.
* `knownSol`: True if the solution is given.
* `numVars`: The number of variables in the system
* `sizeu`: The size of the initial condition (and thus `u`)

### Constructors

SDEProblem(f,σ,u₀;sol=nothing) : Defines the SDE with the specified functions and
defines the solution if sol is given.

"""
type SDEProblem <: DEProblem
  f::Function
  σ::Function
  u₀#::AbstractArray
  sol::Function
  knownSol::Bool
  numVars::Int
  sizeu#::Tuple
  function SDEProblem(f,σ,u₀;sol=nothing)
    if sol==nothing
      knownSol = false
      sol=(u,t,W)->0
    else
      knownSol = true
    end
    if typeof(u₀) <: Number
      sizeu = (1,)
      numVars = 1
    else
      sizeu = size(u₀)
      numVars = size(u₀)[end]
    end
    new(f,σ,u₀,sol,knownSol,numVars,sizeu)
  end
end

"""
ODEProblem

Wraps the data which defines an SDE problem

```math
uₜ = f(uₜ,t)dt
```

with initial condition u₀.

### Fields

* `f`: The drift function in the ODE.
* `u₀`: The initial condition.
* `sol`: A function which describes the solution.
* `knownSol`: True if the solution is given.
* `numVars`: The number of variables in the system
* `sizeu`: The size of the initial condition (and thus `u`)

### Constructors

ODEProblem(f,u₀;sol=nothing) : Defines the SDE with the specified functions and
defines the solution if sol is given.

"""
type ODEProblem <: DEProblem
  f::Function
  u₀#::AbstractArray
  sol::Function
  knownSol::Bool
  numVars::Int
  sizeu#::Tuple
  function ODEProblem(f,u₀;sol=nothing)
    if sol==nothing
      knownSol = false
      sol=(u,t)->0
    else
      knownSol = true
    end
    if typeof(u₀) <: Number
      sizeu = (1,)
      numVars = 1
    else
      sizeu = size(u₀)
      numVars = size(u₀)[end]
    end
    new(f,u₀,sol,knownSol,numVars,sizeu)
  end
end

"""
StokesProblem

Defines the solution to a stationary Stokes problem:

```math

```

### Fields

* `f₁::Function`
* `f₂::Function`
* `g::Function`
* `ugD::Function`
* `vgD::Function`
* `usol::Function`
* `vsol::Function`
* `psol::Function`
* `trueKnown::Bool`

### Constructors

`StokesProblem(f₁,f₂,g,usol,vsol,psol)`

`StokesProblem(f₁,f₂,g,ugD,vgD)`
"""
type StokesProblem
  f₁::Function
  f₂::Function
  g::Function
  ugD::Function
  vgD::Function
  usol::Function
  vsol::Function
  psol::Function
  trueKnown::Bool
  StokesProblem(f₁,f₂,g,usol,vsol,psol) = new(f₁,f₂,g,usol,vsol,usol,vsol,psol,true)
  StokesProblem(f₁,f₂,g,ugD,vgD) = new(f₁,f₂,g,ugD,vgD,nothing,nothing,nothing,false)
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
