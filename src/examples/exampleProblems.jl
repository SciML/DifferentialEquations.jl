### SDE Examples

function linearSDEExample(;α=1,β=1,u₀=1/2)
  f(u,t) = α*u
  σ(u,t) = β*u
  sol(u₀,t,W) = u₀*exp((α-(β^2)/2)*t+β*W)
  return(SDEProblem(f,σ,u₀,sol=sol))
end

function twoDimlinearSDEExample(;α=ones(4,2),β=ones(4,2),u₀=ones(4,2)/2)
  f(u,t) = α.*u
  σ(u,t) = β.*u
  sol(u₀,t,W) = u₀.*exp((α-(β.^2)./2).*t+β.*W)
  return(SDEProblem(f,σ,u₀,sol=sol))
end

function cubicSDEExample(;u₀=1/2)
  f(u,t) = -.25*u.*(1-u.^2)
  σ(u,t) = .5*(1-u.^2)
  sol(u₀,t,W) = ((1+u₀)*exp(W)+u₀-1)./((1+u₀)*exp(W)+1-u₀)
  return(SDEProblem(f,σ,u₀,sol=sol))
end

function waveSDEExample(;u₀=1)
  f(u,t) = -0.01*sin(u).*cos(u).^3
  σ(u,t) = 0.1*cos(u).^2
  sol(u₀,t,W) = atan(0.1*W + tan(u₀))
  return(SDEProblem(f,σ,u₀,sol=sol))
end

function additiveSDEExample(;α=0.1,β=0.5,u₀=1)
  f(u,t) = β./sqrt(1+t) - u./(2*(1+t))
  σ(u,t) = α*β./sqrt(1+t)
  sol(u₀,t,W) = u₀./sqrt(1+t) + β*(t+α*W)./sqrt(1+t)
  return(SDEProblem(f,σ,u₀,sol=sol))
end

function multiDimAdditiveSDEExample(;α=[0.1;0.1;0.1;0.1],β=0.5,u₀=1)
  f(u,t) = β/sqrt(1+t) - u/(2*(1+t))
  σ(u,t) = α*β/sqrt(1+t)
  sol(u₀,t,W) = u₀/sqrt(1+t) + β*(t+α*W)/sqrt(1+t)
  return(SDEProblem(f,σ,u₀,sol=sol))
end
### Finite Element Examples

"Example problem with solution: ``u(x,y,t)=0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x-t+0.5).^2 + (y-t+0.5).^2))``"
function heatProblemExample_moving()
  sol(x,t) = 0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x[:,1]-t+0.5).^2 + (x[:,2]-t+0.5).^2))
  Du(x,t) = -50[sol(x,t).*(0.5-t+x[:,1]) sol(x,t).*(0.5-t+x[:,2])]
  f(x,t) = (-5).*exp((-25).*((3/2)+6.*t.^2+x[:,1]+x[:,1].^2+x[:,2]+x[:,2].^2+(-2).*t.*(3+x[:,1]+
    x[:,2]))).*((-20)+(-100).*t.^2+(-49).*x[:,1]+(-50).*x[:,1].^2+(-49).*x[:,2]+(-50).*
    x[:,2].^2+2.*t.*(47+50.*x[:,1]+50.*x[:,2])+exp(25.*(1+(-2).*t).^2).*(22+
    100.*t.^2+49.*x[:,1]+50.*x[:,1].^2+49.*x[:,2]+50.*x[:,2].^2+(-2).*t.*(49+50.*x[:,1]+50.*x[:,2])))
  return(HeatProblem(sol,Du,f))
end

"Example problem with solution: ``u(x,y,t)=exp(-10((x-.5).^2 + (y-.5).^2 )-t)``"
function heatProblemExample_diffuse()
  sol(x,t) = exp(-10((x[:,1]-.5).^2 + (x[:,2]-.5).^2 )-t)
  f(x,t)   = exp(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
  Du(x,t) = -20[sol(x,t).*(x[:,1]-.5) sol(x,t).*(x[:,2]-.5)]
  return(HeatProblem(sol,Du,f))
end

"Example problem which starts with 1 at (0.5,0.5) and solves with ``f=gD=0``"
function heatProblemExample_pure()
  f(x,t)  = zeros(size(x,1))
  u₀(x) = float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6)) #Only mass at middle of (0,1)^2
  return(HeatProblem(u₀,f))
end

"Example problem which starts with 0 and solves with ``f(u)=1-u/2``"
function heatProblemExample_birthdeath()
  f(u,x,t)  = ones(size(x,1)) - .5u
  u₀(x) = zeros(size(x,1))
  return(HeatProblem(u₀,f))
end

function heatProblemExample_birthdeathsystem()
  f₁(u,x,t)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x,t)  = ones(size(x,1)) -   u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f))
end

function heatProblemExample_diffusionconstants(;D=[.01 .001],max=1)
  f₁(u,x,t)  = zeros(size(x,1))
  f₂(u,x,t)  = zeros(size(x,1))
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = [max*float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6)) max*float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6))]  # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

function heatProblemExample_birthdeathinteractingsystem()
  f₁(u,x,t)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x,t)  = .5u[:,1] -   u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f))
end

function heatProblemExample_grayscott(;ρ=.03,k=.062,D=[1e-3 .5e-3])
  f₁(u,x,t)  = + u[:,1].*u[:,2].*u[:,2] + ρ*(1-u[:,2])
  f₂(u,x,t)  = u[:,1].*u[:,2].*u[:,2] -(ρ+k).*u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = [ones(size(x,1))+rand(size(x,1)) .25.*float(((.2.<x[:,1].<.6) &
          (.2.<x[:,2].<.6)) | ((.85.<x[:,1]) & (.85.<x[:,2])))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

function heatProblemExample_grayscott(;ρ=.03,k=.062,D=[1e-3 .5e-3])
  f₁(u,x,t)  = - u[:,1].*u[:,2].*u[:,2] + ρ*(1-u[:,1])
  f₂(u,x,t)  = u[:,1].*u[:,2].*u[:,2] -(ρ+k).*u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = [ones(size(x,1))+rand(size(x,1)) .25.*float(((.2.<x[:,1].<.6) &
          (.2.<x[:,2].<.6)) | ((.85.<x[:,1]) & (.85.<x[:,2])))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

function heatProblemExample_gierermeinhardt(;a=1,α=1,D=[0.01 1.0],ubar=1,vbar=0,β=10)
  f₁(u,x,t)  = a*u[:,1].*u[:,1]./u[:,2] + ubar - α*u[:,1]
  f₂(u,x,t)  = a*u[:,1].*u[:,1] + vbar -β.*u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  uss = (ubar +β)/α
  vss = (α/β)*uss.^2
  u₀(x) = [uss*ones(size(x,1))+rand(size(x,1)) vss*ones(size(x,1))+rand(size(x,1))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

"Example problem which starts with 0 and solves with ``f(u)=1-u/2`` with noise ``σ(u)=10u^2``"
function heatProblemExample_stochasticbirthdeath()
  f(u,x,t)  = ones(size(x,1)) - .5u
  u₀(x) = zeros(size(x,1))
  σ(u,x,t) = 10u.^2
  return(HeatProblem(u₀,f,σ=σ))
end

"Example problem with solution: ``u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)``"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  return(PoissonProblem(f,sol,Du))
end

"Example problem with deterministic solution: ``u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)``"
function poissonProblemExample_noisyWave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  σ(x) = 5 #Additive noise
  return(PoissonProblem(f,sol,Du,σ=σ))
end

"Example problem for nonlinear Poisson equation. Uses ``f(u)=1-u/2``."
function poissonProblemExample_birthdeath()
  f(u,x)  = ones(size(x,1)) - .5u
  numVars = 1
  return(PoissonProblem(f,numVars=numVars))
end

function poissonProblemExample_birthdeathsystem()
  f₁(u,x)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x)  = ones(size(x,1)) -   u[:,2]
  f(u,x) = [f₁(u,x) f₂(u,x)]
  u₀(x) = .5*ones(size(x,1),2) # size (x,2), 2 meaning 2 variables
  return(PoissonProblem(f,u₀=u₀))
end

function poissonProblemExample_birthdeathinteractingsystem()
  f₁(u,x)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x)  = .5u[:,1] -   u[:,2]
  f(u,x) = [f₁(u,x) f₂(u,x)]
  u₀(x) = ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
  return(PoissonProblem(f,u₀=u₀))
end

## Stokes Examples

function homogeneousStokesExample(;C=0)
  f₁(x,y)   = zeros(x)
  f₂(x,y)   = zeros(x)
  usol(x,y) = 20x.*y.^3
  vsol(x,y) = 5x.^4 - 5y.^4
  psol(x,y) = 60x.^2 .*y - 20y.^3 + C
  g(x,y)    = 0
  return(StokesProblem(f₁,f₂,g,usol,vsol,psol))
end

function dirichletzeroStokesExample(;C=0)
  f₁(x,y)   = zeros(x)
  f₂(x,y)   = zeros(x)
  usol(x,y) = 0
  vsol(x,y) = 0
  psol(x,y) = 0
  g(x,y)    = 0
  return(StokesProblem(f₁,f₂,g,usol,vsol,psol))
end
