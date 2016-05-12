@doc """
`HeatProblem`

Wraps the data that define a 2D linear heat equation problem:

``u_t = Δu + f(x,t)``

#Constructors

HeatProblem(sol,Du,f): Defines the Dirichlet problem with solution sol, solution gradient Du = [u_x,u_y], and f.

HeatProblem(u0,f,gD,gN): Defines the problem with initial value u0 (as a function or vector), f,
Dirichlet boundary function gD, and Neumann boundary function gN.
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
  function HeatProblem(sol,Du,f)
    u0(x) = sol(x,0)
    gD = sol
    return(new(u0,Du,f,gD,nothing,sol,true))
  end
  function HeatProblem(u0,f,gD,gN)
    knownSol = false
    return(new(u0,nothing,f,gD,gN,nothing,false))
  end
end
"""
PoissonProblem

Wraps the data that define a 2D linear Poisson equation problem:

Δu = f(x,t)

#Constructors

PoissonProblem(f,sol,Du,gN): Defines the Dirichlet problem with solution sol, solution gradient Du = [u_x,u_y],
f, and Neumann boundary data gN.

PoissonProblem(u0,f,gD,gN): Defines the problem with initial value u0 (as a function or vector), f,
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
  function PoissonProblem(f,sol,Du,gN)
    return(new(f,sol,Du,sol,gN,true))
  end
  function PoissonProblem(f,gD,gN)
    return(new(f,nothing,nothing,gD,gN,false))
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
  return(HeatProblem(sol,Du,f))
end

"Example problem with solution: u(x,y,t)=exp(-10((x-.5).^2 + (y-.5).^2 )-t)"
function heatProblemExample_diffuse()
  sol(x,t) = exp(-10((x[:,1]-.5).^2 + (x[:,2]-.5).^2 )-t)
  f(x,t)   = exp(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
  Du(x,t) = -20[sol(x,t).*(x[:,1]-.5) sol(x,t).*(x[:,2]-.5)]
  return(HeatProblem(sol,Du,f))
end

"Example problem which starts with 1 at (0.5,0.5) and solves with f=gD=0"
function heatProblemExample_pure()
  gD(x,t) = zeros(size(x,1))
  f(x,t)  = zeros(size(x,1))
  u0(x) = float(abs(x[:,1]-.5 .< 1e-6) & abs(x[:,2]-.5 .< 1e-6)) #Only mass at middle of (0,1)^2
  return(HeatProblem(u0,f,gD,nothing))
end

"Example problem with solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  gN(x) = 0
  return(PoissonProblem(f,sol,Du,gN))
end
