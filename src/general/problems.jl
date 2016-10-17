"""
`HeatProblem`

Wraps the data that define a 2D heat equation problem:

```math
u_t = Δu + f
```

with bounday conditions `gD` on the dirichlet boundary and gN on the neumann boundary.
Linearity is determined by whether the forcing function `f` is a function of two
variables `(t,x)` or three `(t,x,u)` (with `x=[:,1]` and `y=[:,2]`).

If they keyword `σ` is given, then this wraps the data that define a 2D stochastic heat equation

```math
u_t = Δu + f + σdW_t
```

### Constructors

* `HeatProblem(analytic,Du,f)`: Defines the dirichlet problem with solution `analytic`,
  solution gradient `Du = [u_x,u_y]`, and the forcing function `f`.

* `HeatProblem(u₀,f)`: Defines the problem with initial value `u₀` (as a function) and `f`.
  If your initial data is a vector, wrap it as `u₀(x) = vector`.

Note: If all functions are of `(t,x)`, then the program assumes it's linear. Write
your functions using the math to program syntrax translation: ``x`` `= x[:,1]` and ``y`` `= x[:,2]`.
Use `f=f(t,x,u)` and `σ=σ(t,x,u)` (if specified) for nonlinear problems
(with the boundary conditions still (t,x)). Systems of equations can be specified
with `u_i = u[:,i]` as the ith variable. See the example problems for more help.

### Keyword Arguments

* `gD` = dirichlet boundary function

* `gN` = neumann boundary function

* `σ` = The function which multiplies the noise dW. By default `σ=0`.

* `noisetype` = A string which specifies the type of noise to be generated. By default
  `noisetype=:White` for Gaussian Spacetime White Noise.

* `numvars` = Number of variables in the system. Automatically calculated from u₀ in most cases.

* `D` = Array which defines the diffusion coefficients. Default is `D=ones(1,numvars)`.
"""
type HeatProblem <: DEProblem
  "u₀: Initial value function"
  u₀#::Function
  "Du: Function for the solution gradient [u_x,u_y]"
  Du::Function
  "f: Forcing function in heat equation"
  f::Function
  "gD: dirichlet boundary data"
  gD#::Function
  "gN: neumann boundary data"
  gN#::Function
  "analytic: Solution to the heat problem"
  analytic::Function
  "knownanalytic: Boolean which states whether the solution function is given"
  knownanalytic::Bool
  "islinear: Boolean which states whether the problem is linear or nonlinear"
  islinear::Bool
  numvars::Int
  σ::Function
  stochastic::Bool
  noisetype::Symbol
  D#AbstractArray
  function HeatProblem(analytic,Du,f;gN=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
    islinear = numparameters(f)==2
    knownanalytic = true
    u₀(x) = analytic(0,x)
    numvars = size(u₀([0 0
                       0 0
                       0 0]),2)
    gD = analytic
    if gN == nothing
      gN=(t,x)->zeros(size(x,1),numvars)
    end
    if σ==nothing
      stochastic=false
      σ=(t,x)->zeros(size(x,1),numvars)
    else
      stochastic=true
    end
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
    return(new(u₀,Du,f,gD,gN,analytic,knownanalytic,islinear,numvars,σ,stochastic,noisetype,D))
  end
  function HeatProblem(u₀,f;gD=nothing,gN=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
    if σ==nothing
      stochastic=false
      σ=(t,x)->zeros(size(x,1))
    else
      stochastic=true
    end
    islinear = numparameters(f)==2
    knownanalytic = false
    if islinear
      if u₀==nothing
        u₀=(x)->zeros(size(x,1))
      end
      if gD == nothing
        gD=(t,x)->zeros(size(x,1))
      end
      if gN == nothing
        gN=(t,x)->zeros(size(x,1))
      end
      if D == nothing
        D = 1.0
      end
      numvars = 1
    end
    if !islinear #nonlinear
      if u₀==nothing && numvars == nothing
        warn("u₀ and numvars must be given. numvars assumed 1.")
        numvars = 1
        u₀=(x)->zeros(size(x,1),numvars)
        if gD == nothing
          gD=(t,x)->zeros(size(x,1),numvars)
        end
        if gN == nothing
          gN=(t,x)->zeros(size(x,1),numvars)
        end
        if D == nothing
          D = 1.0
        end
      elseif u₀==nothing #numvars!=nothing
        u₀=(x)->zeros(size(x,1),numvars) #Default to zero
        if gD == nothing
          gD=(t,x)->zeros(size(x,1),numvars)
        end
        if gN == nothing
          gN=(t,x)->zeros(size(x,1),numvars)
        end
        if D == nothing
          D = ones(1,numvars)
        end
      elseif numvars==nothing #If u₀ is given but numvars is not, we're still okay. Generate from size in function.
        numvars=0 #Placeholder, update gD and gN in solver
      end
    end
    return(new(u₀,(x)->0,f,gD,gN,(x)->0,knownanalytic,islinear,numvars,σ,stochastic,noisetype,D))
  end
end

doc"""
`PoissonProblem`

Wraps the data that define a 2D linear Poisson equation problem:

```math
-Δu = f
```

with bounday conditions `gD` on the dirichlet boundary and gN on the neumann boundary.
Linearity is determined by whether the forcing function `f` is a function of one
variable `(x)` or two `(u,x)` (with `x=[:,1]` and `y=[:,2]`).

If they keyword `σ` is given, then this wraps the data that define a 2D stochastic heat equation

```math
-Δu = f + σdW
```

### Constructors

`PoissonProblem(f,analytic,Du)`: Defines the dirichlet problem with analytical solution `analytic`, solution gradient `Du = [u_x,u_y]`,
and forcing function `f`

`PoissonProblem(u₀,f)`: Defines the problem with initial value `u₀` (as a function) and f.
If your initial data is a vector, wrap it as `u₀(x) = vector`.

Note: If all functions are of `(x)`, then the program assumes it's linear. Write
your functions using the math to program syntrax translation: ``x`` `= x[:,1]` and ``y`` `= x[:,2]`.
Use `f=f(u,x)` and `σ=σ(u,x)` (if specified) for nonlinear problems
(with the boundary conditions still (x)). Systems of equations can be specified
with `u_i = u[:,i]` as the ith variable. See the example problems for more help.

### Keyword Arguments

* `gD` = dirichlet boundary function

* `gN` = neumann boundary function

* `σ` = The function which multiplies the noise ``dW``. By default `σ=0`.

* `noisetype` = A string which specifies the type of noise to be generated. By default
  `noisetype=:White` for Gaussian Spacetime White Noise.

* `numvars` = The number of variables in the Poisson system. Automatically calculated in many cases.

* `D` = Vector of diffusion coefficients. Defaults is `D=ones(1,numvars)`.

"""
type PoissonProblem <: DEProblem
  "f: Forcing function in the Poisson problem"
  f#::Function
  "analytic: Solution to the Poisson problem"
  analytic::Function
  "Du: Gradient of the solution to the Poisson problem"
  Du::Function
  "gD: dirichlet Boundary Data"
  gD#::Nullable{Function}
  "gN: neumann Boundary Data"
  gN#::Nullable{Function}
  "knownanalytic: Boolean which states whether the solution function is given"
  knownanalytic::Bool
  "islinear: Boolean which states whether the problem is linear or nonlinear"
  islinear::Bool
  u₀::Function
  numvars::Int
  σ::Function
  stochastic::Bool
  noisetype::Symbol
  D#::AbstractArray
  function PoissonProblem(f,analytic,Du;gN=nothing,σ=nothing,u₀=nothing,noisetype=:White,numvars=nothing,D=nothing)
    gD = analytic
    numvars = size(analytic([0 0
                        0 0
                        0 0]),2)
    islinear = numparameters(f)==1
    if gN == nothing
      gN=(x)->zeros(size(x,1),numvars)
    end
    if u₀==nothing
      u₀=(x)->zeros(size(x,1),numvars)
    end
    if D == nothing
      if numvars == 1
        D = 1.0
      else
        D = ones(1,numvars)
      end
    end
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1),numvars)
    else
      stochastic=true
    end
    return(new(f,analytic,Du,analytic,gN,true,islinear,u₀,numvars,σ,stochastic,noisetype,D))
  end
  function PoissonProblem(f;gD=nothing,gN=nothing,u₀=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
    if σ==nothing
      stochastic=false
      σ=(x)->zeros(size(x,1))
    else
      stochastic = true
    end
    islinear = numparameters(f)==1
    if islinear && u₀==nothing
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
      numvars = 1
    end
    if !islinear #nonlinear
      if u₀==nothing && numvars == nothing
        warn("u₀ and numvars must be given. numvars assumed 1.")
        numvars = 1
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
      elseif u₀==nothing #numvars!=nothing
        u₀=(x)->zeros(size(x,1),numvars) #Default to zero
        if gD == nothing
          gD=(x)->zeros(size(x,1),numvars)
        end
        if gN == nothing
          gN=(x)->zeros(size(x,1),numvars)
        end
        if D == nothing
          D = ones(1,numvars)
        end
      elseif numvars==nothing #If u₀ is given but numvars is not, we're still okay. Generate from size in function.
        numvars=0 #Placeholder, update gD and gN in solver
      end
    end
    return(new(f,(x)->0,(x)->0,gD,gN,false,islinear,u₀,numvars,σ,stochastic,noisetype,D))
  end
end

"""
`SDEProblem`

Wraps the data which defines an SDE problem

```math
u = f(u,t)dt + Σgᵢ(u,t)dWⁱ
```

with initial condition ``u₀``.

### Constructors

`SDEProblem(f,g,u₀;analytic=nothing)` : Defines the SDE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the SDE.
* `g`: The noise function in the SDE.
* `u₀`: The initial condition.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system
* `sizeu`: The size of the initial condition (and thus `u`)
* `noise`: The noise process applied to the noise upon generation.

"""
type SDEProblem <: AbstractSDEProblem
  f::Function
  g::Function
  u₀#::AbstractArray
  analytic::Function
  knownanalytic::Bool
  numvars::Int
  sizeu#::Tuple
  isinplace::Bool
  noise::NoiseProcess
  function SDEProblem(f,g,u₀;analytic=nothing,noise=WHITE_NOISE)
    isinplace = numparameters(f)>=3
    if analytic==nothing
      knownanalytic = false
      analytic=(t,u,W)->0
    else
      knownanalytic = true
    end
    if typeof(u₀) <: Number
      sizeu = (1,)
      numvars = 1
    else
      sizeu = size(u₀)
      numvars = size(u₀)[end]
    end
    new(f,g,u₀,analytic,knownanalytic,numvars,sizeu,isinplace,noise)
  end
end

"""
`ODEProblem`

Wraps the data which defines an ODE problem

```math
\\frac{du}{dt} = f(t,u)
```

with initial condition ``u₀``.

### Constructors

`ODEProblem(f,u₀;analytic=nothing)` : Defines the ODE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the ODE.
* `u₀`: The initial condition.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system

"""
type ODEProblem{uType<:Union{AbstractArray,Number},uEltype<:Number} <: AbstractODEProblem
  f::Function
  u₀::uType
  analytic::Function
  knownanalytic::Bool
  numvars::Int
  isinplace::Bool
end

function ODEProblem(f::Function,u₀;analytic=nothing)
  isinplace = numparameters(f)>=3
  if analytic==nothing
    knownanalytic = false
    analytic=(t,u,du)->0
  else
    knownanalytic = true
  end
  if typeof(u₀) <: Number
    sizeu = (1,)
    numvars = 1
  else
    numvars = size(u₀)[end]
  end
  ODEProblem{typeof(u₀),eltype(u₀)}(f,u₀,analytic,knownanalytic,numvars,isinplace)
end

"""
`DAEProblem`

Wraps the data which defines an DAE problem

```math
f(t,u,du₀) = 0
```

with initial conditions ``u₀`` and ``du₀``.

### Constructors

`DAEProblem(f,u₀,du₀;analytic=nothing)` : Defines the SDE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the DAE.
* `u₀`: The initial condition.
* `du₀`: The initial derivative.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system

"""
type DAEProblem{uType<:Union{AbstractArray,Number},uEltype<:Number,rateType<:Union{AbstractArray,Number}} <: AbstractODEProblem
  f::Function
  u₀::uType
  du₀::rateType
  analytic::Function
  knownanalytic::Bool
  numvars::Int
  isinplace::Bool
end

function DAEProblem(f::Function,u₀,du₀;analytic=nothing)
  isinplace = numparameters(f)>=4
  if analytic==nothing
    knownanalytic = false
    analytic=(t,u,du,out)->0
  else
    knownanalytic = true
  end
  if typeof(u₀) <: Number
    sizeu = (1,)
    numvars = 1
  else
    numvars = size(u₀)[end]
  end
  DAEProblem{typeof(u₀),eltype(u₀),typeof(du₀)}(f,u₀,du₀,analytic,knownanalytic,numvars,isinplace)
end


"""
`StokesProblem`

Defines the solution to a stationary Stokes problem:

```math

```

### Constructors

`StokesProblem(f₁,f₂,g,uanalytic,vanalytic,panalytic)`

`StokesProblem(f₁,f₂,g,ugD,vgD)`

### Fields

* `f₁::Function`
* `f₂::Function`
* `g::Function`
* `ugD::Function`
* `vgD::Function`
* `uanalytic::Function`
* `vanalytic::Function`
* `panalytic::Function`
* `trueknown::Bool`
"""
type StokesProblem
  f₁#::Function
  f₂#::Function
  g#::Function
  ugD#::Function
  vgD#::Function
  uanalytic#::Function
  vanalytic#::Function
  panalytic#::Function
  trueknown::Bool
  StokesProblem(f₁,f₂,g,uanalytic,vanalytic,panalytic) = new(f₁,f₂,g,uanalytic,vanalytic,uanalytic,vanalytic,panalytic,true)
  StokesProblem(f₁,f₂,g,ugD,vgD) = new(f₁,f₂,g,ugD,vgD,nothing,nothing,nothing,false)
end
