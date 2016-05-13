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

#=
@doc "Example problem with deterministic solution: u(x,y,t)=0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x-t+0.5).^2 + (y-t+0.5).^2))" ->
function heatProblemExample_stochasticMoving()
  sol(x,t) = 0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x[:,1]-t+0.5).^2 + (x[:,2]-t+0.5).^2))
  Du(x,t) = -50[sol(x,t).*(0.5-t+x[:,1]) sol(x,t).*(0.5-t+x[:,2])]
  f(u,x,t) = (-5).*exp((-25).*((3/2)+6.*t.^2+x[:,1]+x[:,1].^2+x[:,2]+x[:,2].^2+(-2).*t.*(3+x[:,1]+
    x[:,2]))).*((-20)+(-100).*t.^2+(-49).*x[:,1]+(-50).*x[:,1].^2+(-49).*x[:,2]+(-50).*
    x[:,2].^2+2.*t.*(47+50.*x[:,1]+50.*x[:,2])+exp(25.*(1+(-2).*t).^2).*(22+
    100.*t.^2+49.*x[:,1]+50.*x[:,1].^2+49.*x[:,2]+50.*x[:,2].^2+(-2).*t.*(49+50.*x[:,1]+50.*x[:,2])))
  isLinear = false
  stochastic = true
  σ(u,x,t) = 1000u.*x[:,1].*(1-x[:,1]).*x[:,2].*(1-x[:,2])
  return(HeatProblem(sol,Du,f,isLinear,σ=σ,stochastic=stochastic))
end
# Bad example, too stiff on f to see noise...
=#

"Example problem with solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  gN(x) = 0
  isLinear = true
  return(PoissonProblem(f,sol,Du,gN,isLinear))
end

"Example problem with deterministic solution: u(x,y,t)= sin(2π.*x).*cos(2π.*y)/(8π*π)"
function poissonProblemExample_noisyWave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  gN(x) = 0
  isLinear = true
  stochastic = true
  σ(x) = 5 #Additive noise
  return(PoissonProblem(f,sol,Du,gN,isLinear,σ=σ,stochastic=stochastic))
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
