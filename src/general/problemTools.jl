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
  function HeatProblem(sol,Du,f,isLinear)
    knownSol = true
    u0(x) = sol(x,0)
    gD = sol
    return(new(u0,Du,f,gD,nothing,sol,knownSol,isLinear))
  end
  function HeatProblem(u0,f,gD,gN,isLinear)
    knownSol = false
    return(new(u0,nothing,f,gD,gN,nothing,knownSol,isLinear))
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
  function PoissonProblem(f,sol,Du,gN,isLinear)
    return(new(f,sol,Du,sol,gN,true,isLinear))
  end
  function PoissonProblem(f,gD,gN,isLinear)
    return(new(f,nothing,nothing,gD,gN,false,isLinear))
  end
end

@doc "Example problem with solution: u(x,y,t)=0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x-t+0.5).^2 + (y-t+0.5).^2))" ->
function heatProblemExample_moving()
  sol(x,t) = 0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x[:,1]-t+0.5).^2 + (x[:,2]-t+0.5).^2))
  Du(x,t) = -50[sol(x,t).*(0.5-t+x[:,1]) sol(x,t).*(0.5-t+x[:,2])]
  f(x,t) = (-5).*exp((-25).*((3/2)+6.*t.^2+x[:,1]+x[:,1].^2+x[:,2]+x[:,2].^2+(-2).*t.*(3+x[:,1]+
    x[:,2]))).*((-20)+(-100).*t.^2+(-49).*x[:,1]+(-50).*x[:,1].^2+(-49).*x[:,2]+(-50).*
    x[:,2].^2+2.*t.*(47+50.*x[:,1]+50.*x[:,2])+exp(25.*(1+(-2).*t).^2).*(22+
    100.*t.^2+49.*x[:,1]+50.*x[:,1].^2+49.*x[:,2]+50.*x[:,2].^2+(-2).*t.*(49+50.*x[:,1]+50.*x[:,2])))
  isLinear = true
  return(HeatProblem(sol,Du,f,isLinear))
end

"Example problem with solution: u(x,y,t)=exp(-10((x-.5).^2 + (y-.5).^2 )-t)"
function heatProblemExample_diffuse()
  sol(x,t) = exp(-10((x[:,1]-.5).^2 + (x[:,2]-.5).^2 )-t)
  f(x,t)   = exp(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
  Du(x,t) = -20[sol(x,t).*(x[:,1]-.5) sol(x,t).*(x[:,2]-.5)]
  isLinear = true
  return(HeatProblem(sol,Du,f,isLinear))
end

"Example problem which starts with 1 at (0.5,0.5) and solves with f=gD=0"
function heatProblemExample_pure()
  gD(x,t) = zeros(size(x,1))
  f(x,t)  = zeros(size(x,1))
  u0(x) = float(abs(x[:,1]-.5 .< 1e-6) & abs(x[:,2]-.5 .< 1e-6)) #Only mass at middle of (0,1)^2
  isLinear = true
  return(HeatProblem(u0,f,gD,isLinear))
end

"Example problem which starts with 0 and solves with f(u)=1-.1u"
function heatProblemExample_birthdeath()
  gD(x,t) = zeros(size(x,1))
  f(u,x,t)  = ones(size(x,1)) - .5u
  u0(x) = zeros(size(x,1))
  gN(x,t) = 0
  isLinear = false
  return(HeatProblem(u0,f,gD,gN,isLinear))
end

"Example problem with solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  gN(x) = 0
  isLinear = true
  return(PoissonProblem(f,sol,Du,gN,isLinear))
end

function poissonProblemExample_birthdeath()
  gD(x) = 0
  f(u,x)  = ones(size(x,1)) - .5u
  gN(x) = 0
  isLinear = false
  return(PoissonProblem(f,gD,gN,isLinear))
end

#=
function poissonProblemExample_nonlinearPBE()
  f(u,x) = -sinh(u)*(h^2)
  ubar(s) = log((1+cos(s))/(1-cos(s)))
  function sol(x)
    N = size(x)[1]
    res = Array{Float64}(N,N) #Assumes square
    a = [1.,2.]/sqrt(5)
    for i = 1:N, j = 1:N #Assumes square
      z = 0.1 + x[i,1]*a[1] + x[j,2]*a[2]
      res[i,j] = ubar(z)
    end
    return(res)
  end
  Du(x) =
end
=#
