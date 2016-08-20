### ODE Examples

"""Example problem with solution ``u(t)=u₀*exp(α*t)``"""
function linearODEExample(;α=1,u₀=1/2)
  f(u,t) = α*u
  analytic(u₀,t) = u₀*exp(α*t)
  return(ODEProblem(f,u₀,analytic=analytic))
end

"""Van der Pol Equations. For difficult version, use α=1e6"""
function vanDerPolExample(α=1,u₀=[0;sqrt(3)])
  function f(du,u,t)
    du[1] = ((1-u[2].^2)*u[1] - u[2])*α
    du[2] = u[1]
  end
  return(ODEProblem(f,u₀))
end

"""Example problem of 8 linear ODEs (as a 4x2 matrix) with solution ``u(t)=exp(α.*t)`` and random initial conditions"""
function twoDimlinearODEExample(;α=ones(4,2),u₀=rand(4,2).*ones(4,2)/2)
  f(u,t) = α.*u
  analytic(u₀,t) = u₀.*exp(α.*t)
  return(ODEProblem(f,u₀,analytic=analytic))
end

"""Example problem of 8 linear ODEs (as a 4x2 matrix) with solution ``u(t)=exp(α.*t)`` and random initial conditions"""
function twoDimlinearODEExample!(;α=ones(4,2),u₀=rand(4,2).*ones(4,2)/2)
  function f(du,u,t)
    @inbounds for i in eachindex(u)
      du[i] = α[i]*u[i]
    end
    du
  end
  analytic(u₀,t) = u₀.*exp(α.*t)
  return(ODEProblem(f,u₀,analytic=analytic))
end

function lorenzAttractorODEExample(;σ=10.,ρ=28.,β=8//3,u₀=ones(3))
  f₁(u,t) = σ*(u[2]-u[1])
  f₂(u,t) = u[1]*(ρ-u[3]) - u[2]
  f₃(u,t) = u[1]*u[2] - β*u[3]
  f(u,t) = [f₁(u,t);f₂(u,t);f₃(u,t)]
  return(ODEProblem(f,u₀))
end

function lorenzAttractorODEExample!(;σ=10.,ρ=28.,β=8//3,u₀=ones(3))
  function f(du,u,t)
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
  end
  return(ODEProblem(f,u₀))
end

"""
The Robertson biochemical reactions

http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Single/DemoRobertson/demorobertson.pdf

Or from Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 129

Usually solved on [0,1e11]
"""
function ROBERODEExample(u₀=[1.0;0.0;0.0],k₁=0.04,k₂=3e7,k₃=1e4)
  function f(du,u,t)
    du[1] = -k₁*u[1]+k₃*u[2]*u[3]
    du[2] =  k₁*u[1]-k₂*(u[2])^2-k₃*u[2]*u[3]
    du[3] =  k₂*(u[2])^2
    nothing
  end
  return(ODEProblem(f,u₀))
end

"""
From Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 129

Usually solved on t₀ = 0.0; T = parse(BigFloat,"17.0652165601579625588917206249")
Periodic with that setup
"""
function threebodyODEExample(u₀=[0.994, 0.0, 0.0, parse(BigFloat,"-2.00158510637908252240537862224")])
  μ = parse(BigFloat,"0.012277471"); μ′ = 1 - μ
  ToT = 3/2
  function f(du,u,t)
    # 1 = y₁
    # 2 = y₂
    # 3 = y₁'
    # 4 = y₂'
    D₁ = ((u[1]+μ)^2 + u[2]^2)^ToT
    D₂ = ((u[1]-μ′)^2 + u[2]^2)^ToT
    du[1] = u[3]
    du[2] = u[4]
    du[3] = u[1] + 2u[4] - μ′*(u[1]+μ)/D₁ - μ*(u[1]-μ′)/D₂
    du[4] = u[2] - 2u[3] - μ′*u[2]/D₁ - μ*u[2]/D₂
    nothing
  end
  return(ODEProblem(f,u₀))
end

### SDE Examples

"""Example problem with solution ``u(t,W)=u₀*exp((α-(β^2)/2)*t+β*W)``"""
function linearSDEExample(;α=1,β=1,u₀=1/2)
  f(u,t) = α*u
  σ(u,t) = β*u
  analytic(u₀,t,W) = u₀*exp((α-(β^2)/2)*t+β*W)
  return(SDEProblem(f,σ,u₀,analytic=analytic))
end

"""Example problem of 8 linear SDEs (as a 4x2 matrix) with solution ``u(t,W)=u₀*exp((α-(β^2)/2)*t+β*W)``"""
function twoDimlinearSDEExample(;α=ones(4,2),β=ones(4,2),u₀=ones(4,2)/2)
  f(u,t) = α.*u
  σ(u,t) = β.*u
  analytic(u₀,t,W) = u₀.*exp((α-(β.^2)./2).*t+β.*W)
  return(SDEProblem(f,σ,u₀,analytic=analytic))
end

"""Example problem with solution ``u(t,W)=((1+u₀)*exp(W)+u₀-1)./((1+u₀)*exp(W)+1-u₀)``"""
function cubicSDEExample(;u₀=1/2)
  f(u,t) = -.25*u.*(1-u.^2)
  σ(u,t) = .5*(1-u.^2)
  analytic(u₀,t,W) = ((1+u₀)*exp(W)+u₀-1)./((1+u₀)*exp(W)+1-u₀)
  return(SDEProblem(f,σ,u₀,analytic=analytic))
end

"""Example problem with solution ``u(t,W)=atan(0.1*W + tan(u₀))``"""
function waveSDEExample(;u₀=1.)
  f(u,t) = -0.01*sin(u).*cos(u).^3
  σ(u,t) = 0.1*cos(u).^2
  analytic(u₀,t,W) = atan(0.1*W + tan(u₀))
  return(SDEProblem(f,σ,u₀,analytic=analytic))
end

"""Example additive noise problem with solution ``u₀./sqrt(1+t) + β*(t+α*W)./sqrt(1+t)``"""
function additiveSDEExample(;α=0.1,β=0.05,u₀=1.)
  f(u,t) = β./sqrt(1+t) - u./(2*(1+t))
  σ(u,t) = α*β./sqrt(1+t)
  analytic(u₀,t,W) = u₀./sqrt(1+t) + β*(t+α*W)./sqrt(1+t)
  return(SDEProblem(f,σ,u₀,analytic=analytic))
end

"""Multiple Ito dimension extension of additiveSDEExample"""
function multiDimAdditiveSDEExample(;α=[0.1;0.1;0.1;0.1],β=[0.5;0.25;0.125;0.1115],u₀=[1.;1.;1.;1.])
  f(u,t) = β./sqrt(1+t) - u./(2*(1+t))
  σ(u,t) = α.*β./sqrt(1+t)
  analytic(u₀,t,W) = u₀./sqrt(1+t) + β.*(t+α.*W)/sqrt(1+t)
  return(SDEProblem(f,σ,u₀,analytic=analytic))
end

function lorenzAttractorSDEExample(;α=10.,ρ=28.,β=8//3,u₀=ones(3),σ₀=1)
  f₁(u,t) = α*(u[2]-u[1])
  f₂(u,t) = u[1]*(ρ-u[3]) - u[2]
  f₃(u,t) = u[1]*u[2] - β*u[3]
  f(u,t) = [f₁(u,t);f₂(u,t);f₃(u,t)]
  σ(u,t) = σ₀ #Additive
  return(SDEProblem(f,σ,u₀))
end

function oval2ModelExample(;largeFluctuations=false,useBigs=false,noiseLevel=1)
  #Parameters
  J1_200=3.
  J1_34=0.15
  J2_200=0.2
  J2_34=0.35
  J_2z=0.9
  J_O=0.918
  J_SO=0.5
  # J_ZEB=0.06
  J_ecad1=0.1
  J_ecad2=0.3
  J_ncad1=0.4
  J_ncad2=0.4
  J_ncad3 = 2
  J_snail0=0.6
  J_snail1=1.8
  J_zeb=3.0
  K1=1.0
  K2=1.0
  K3=1.0
  K4=1.0
  K5=1.0
  KTGF=20.
  Ks=100.
  # TGF0=0
  TGF_flg=0.
  Timescale=1000.
  dk_ZR1=0.5
  dk_ZR2=0.5
  dk_ZR3=0.5
  dk_ZR4=0.5
  dk_ZR5=0.5
  k0O=0.35
  k0_200=0.0002
  k0_34=0.001
  k0_snail=0.0005
  k0_zeb=0.003
  kO=1.2
  kOp=10.
  k_200=0.02
  k_34=0.019
  k_OT=1.1
  k_SNAIL=16.
  k_TGF=1.5
  k_ZEB=16.
  k_ecad0=5.
  k_ecad1=15.
  k_ecad2=5.
  k_ncad0=5.
  k_ncad1=2.
  k_ncad2=5.
  k_snail=0.05
  k_tgf=0.05
  k_zeb=0.06
  kdO=1.
  kd_200=0.035
  kd_34=0.035
  kd_SNAIL=1.6
  kd_SR1=0.9
  kd_TGF=0.9
  kd_ZEB=1.66
  kd_ecad=0.05
  kd_ncad=0.05
  kd_snail=0.09
  kd_tgf=0.1
  kd_tgfR=1.0
  kd_zeb=0.1
  kd_Op = 10.
  lamda1=0.5
  lamda2=0.5
  lamda3=0.5
  lamda4=0.5
  lamda5=0.5
  lamdas=0.5
  lamdatgfR=0.8
  nO=6.
  nSO=2.
  nzo=2.
  GE = 1.
  function f(y,t)
    # y(1) = snailt
    # y(2) = SNAIL
    # y(3) = miR34t
    # y(4) = SR1 # abundance of SNAIL/miR34 complex
    # y(5) = zebt
    # y(6) = ZEB
    # y(7) = miR200t
    # y(8) = ZR1 # abundance of ZEB/miR200 complex with i copies of miR200 bound on the sequence of ZEB1
    # y(9) = ZR2
    # y(10) = ZR3
    # y(11) = ZR4
    # y(12) = ZR5
    # y(13) = tgft
    # y(14) = TGF
    # y(15) = tgfR # abundance of TGF/miR200 complex
    # y(16) = Ecad
    # y(17) = Ncad
    # y(18) = Ovol2
    dy = Vector{Float64}(19)
    TGF0=.5(t>100)
    #ODEs
    dy[1]=k0_snail+k_snail*(((y[14]+TGF0)/J_snail0))^2/(1+(((y[14]+TGF0)/J_snail0))^2+(y[19]/J_SO)^nSO)/(1+y[2]/J_snail1)-kd_snail*(y[1]-y[4])-kd_SR1*y[4]
    dy[2]=k_SNAIL*(y[1]-y[4])-kd_SNAIL*y[2]
    dy[3]=k0_34+k_34/(1+((y[2]/J1_34))^2+((y[6]/J2_34))^2)-kd_34*(y[3]-y[4])-kd_SR1*y[4]+lamdas*kd_SR1*y[4]
    dy[4]=Timescale*(Ks*(y[1]-y[4])*(y[3]-y[4])-y[4])
    dy[5]=k0_zeb+k_zeb*((y[2]/J_zeb))^2/(1+((y[2]/J_zeb))^2+((y[19]/J_2z))^nO)-kd_zeb*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-dk_ZR1*5*y[8]-dk_ZR2*10*y[9]-dk_ZR3*10*y[10]-dk_ZR4*5*y[11]-dk_ZR5*y[12]
    dy[6]=k_ZEB*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-kd_ZEB*y[6]
    dy[7]=k0_200+k_200/(1+((y[2]/J1_200))^3+((y[6]/J2_200))^2)-kd_200*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])-dk_ZR1*5*y[8]-dk_ZR2*2*10*y[9]-dk_ZR3*3*10*y[10]-dk_ZR4*4*5*y[11]-dk_ZR5*5*y[12]+lamda1*dk_ZR1*5*y[8]+lamda2*dk_ZR2*2*10*y[9]+lamda3*dk_ZR3*3*10*y[10]+lamda4*dk_ZR4*4*5*y[11]+lamda5*dk_ZR5*5*y[12]-kd_tgfR*y[15]+lamdatgfR*kd_tgfR*y[15]
    dy[8]=Timescale*(K1*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-y[8])
    dy[9]=Timescale*(K2*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[8]-y[9])
    dy[10]=Timescale*(K3*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[9]-y[10])
    dy[11]=Timescale*(K4*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[10]-y[11])
    dy[12]=Timescale*(K5*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[11]-y[12])
    dy[13]=k_tgf-kd_tgf*(y[13]-y[15])-kd_tgfR*y[15]
    dy[14]=k_OT+k_TGF*(y[13]-y[15])-kd_TGF*y[14]
    dy[15]=Timescale*(TGF_flg+KTGF*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*(y[13]-y[15])-y[15])
    dy[16]=GE*(k_ecad0+k_ecad1/(((y[2]/J_ecad1))^2+1)+k_ecad2/(((y[6]/J_ecad2))^2+1)-kd_ecad*y[16])
    dy[17]=k_ncad0+k_ncad1*(((y[2]/J_ncad1))^2)/(((y[2]/J_ncad1))^2+1)+k_ncad2*(((y[6]/J_ncad2))^2)/(((y[6]/J_ncad2)^2+1)*(1+y[19]/J_ncad3))-kd_ncad*y[17]
    dy[18]=k0O+kO/(1+((y[6]/J_O))^nzo)-kdO*y[18]
    dy[19]=kOp*y[18]-kd_Op*y[19]
    return(dy)
  end

  function σ1(y,t)
    dσ = zeros(19)
    dσ[1] = noiseLevel*1.5y[1]
    dσ[18]= noiseLevel*6y[18]
    return(dσ)
  end

  function σ2(y,t)
    dσ = zeros(19)
    dσ[1] = 0.02y[1]
    dσ[16]= 0.02y[16]
    dσ[18]= 0.2y[18]
    dσ[17]= 0.02y[17]
    return(dσ)
  end

  if largeFluctuations
    σ = σ1
  else
    σ = σ2
  end

  u₀ = [0.128483;1.256853;0.0030203;0.0027977;0.0101511;0.0422942;0.2391346;0.0008014;0.0001464;2.67e-05;4.8e-6;9e-7;0.0619917;1.2444292;0.0486676;199.9383546;137.4267984;1.5180203;1.5180203] #Fig 9B
  if useBigs
    u₀ = big(u₀)
  end
  #u₀ =  [0.1701;1.6758;0.0027;0.0025;0.0141;0.0811;0.1642;0.0009;0.0001;0.0000;0.0000;0.0000;0.0697;1.2586;0.0478;194.2496;140.0758;1.5407;1.5407] #Fig 9A
  return(SDEProblem(f,σ,u₀))
end

### Finite Element Examples

"Example problem with solution: ``u(x,y,t)=0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x-t+0.5).^2 + (y-t+0.5).^2))``"
function heatProblemExample_moving()
  analytic(x,t) = 0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x[:,1]-t+0.5).^2 + (x[:,2]-t+0.5).^2))
  Du(x,t) = -50[analytic(x,t).*(0.5-t+x[:,1]) analytic(x,t).*(0.5-t+x[:,2])]
  f(x,t) = (-5).*exp((-25).*((3/2)+6.*t.^2+x[:,1]+x[:,1].^2+x[:,2]+x[:,2].^2+(-2).*t.*(3+x[:,1]+
    x[:,2]))).*((-20)+(-100).*t.^2+(-49).*x[:,1]+(-50).*x[:,1].^2+(-49).*x[:,2]+(-50).*
    x[:,2].^2+2.*t.*(47+50.*x[:,1]+50.*x[:,2])+exp(25.*(1+(-2).*t).^2).*(22+
    100.*t.^2+49.*x[:,1]+50.*x[:,1].^2+49.*x[:,2]+50.*x[:,2].^2+(-2).*t.*(49+50.*x[:,1]+50.*x[:,2])))
  return(HeatProblem(analytic,Du,f))
end

"Example problem with solution: ``u(x,y,t)=exp(-10((x-.5).^2 + (y-.5).^2 )-t)``"
function heatProblemExample_diffuse()
  analytic(x,t) = exp(-10((x[:,1]-.5).^2 + (x[:,2]-.5).^2 )-t)
  f(x,t)   = exp(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
  Du(x,t) = -20[analytic(x,t).*(x[:,1]-.5) analytic(x,t).*(x[:,2]-.5)]
  return(HeatProblem(analytic,Du,f))
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

"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=1-v``"
function heatProblemExample_birthdeathsystem()
  f₁(u,x,t)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x,t)  = ones(size(x,1)) -   u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f))
end

"Example problem which solves the homogeneous Heat equation with all mass starting at (1/2,1/2) with two different diffusion constants."
function heatProblemExample_diffusionconstants(;D=[.01 .001],max=1)
  f₁(u,x,t)  = zeros(size(x,1))
  f₂(u,x,t)  = zeros(size(x,1))
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = [max*float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6)) max*float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6))]  # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=.5u-v``"
function heatProblemExample_birthdeathinteractingsystem()
  f₁(u,x,t)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x,t)  = .5u[:,1] -   u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f))
end

"Example problem which solves the Gray-Scott equations with quasi-random initial conditions"
function heatProblemExample_grayscott(;ρ=.03,k=.062,D=[1e-3 .5e-3])
  f₁(u,x,t)  = + u[:,1].*u[:,2].*u[:,2] + ρ*(1-u[:,2])
  f₂(u,x,t)  = u[:,1].*u[:,2].*u[:,2] -(ρ+k).*u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  u₀(x) = [ones(size(x,1))+rand(size(x,1)) .25.*float(((.2.<x[:,1].<.6) &
          (.2.<x[:,2].<.6)) | ((.85.<x[:,1]) & (.85.<x[:,2])))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

"Example problem which solves the Gierer-Meinhardt equations wtih quasi-random initial perturbations."
function heatProblemExample_gierermeinhardt(;a=1,α=1,D=[0.01 1.0],ubar=1,vbar=0,β=10,startNoise=0.01)
  f₁(u,x,t)  = a*u[:,1].*u[:,1]./u[:,2] + ubar - α*u[:,1]
  f₂(u,x,t)  = a*u[:,1].*u[:,1] + vbar -β.*u[:,2]
  f(u,x,t) = [f₁(u,x,t) f₂(u,x,t)]
  uss = (ubar +β)/α
  vss = (α/β)*uss.^2
  u₀(x) = [uss*ones(size(x,1))+startNoise*rand(size(x,1)) vss*ones(size(x,1))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

"Example problem which starts with 0 and solves with ``f(u)=1-u/2`` with noise ``σ(u)=10u^2``"
function heatProblemExample_stochasticbirthdeath()
  f(u,x,t)  = ones(size(x,1)) - .5u
  u₀(x) = zeros(size(x,1))
  σ(u,x,t) = 1u.^2
  return(HeatProblem(u₀,f,σ=σ))
end

"Example problem with solution: ``u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)``"
function poissonProblemExample_wave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  analytic(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  return(PoissonProblem(f,analytic,Du))
end

"Example problem with deterministic solution: ``u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)``"
function poissonProblemExample_noisyWave()
  f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
  analytic(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
  Du(x) = [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
  σ(x) = 5 #Additive noise
  return(PoissonProblem(f,analytic,Du,σ=σ))
end

"Example problem for nonlinear Poisson equation. Uses ``f(u)=1-u/2``."
function poissonProblemExample_birthdeath()
  f(u,x)  = ones(size(x,1)) - .5u
  numvars = 1
  return(PoissonProblem(f,numvars=numvars))
end

"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=1-v``"
function poissonProblemExample_birthdeathsystem()
  f₁(u,x)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x)  = ones(size(x,1)) -   u[:,2]
  f(u,x) = [f₁(u,x) f₂(u,x)]
  u₀(x) = .5*ones(size(x,1),2) # size (x,2), 2 meaning 2 variables
  return(PoissonProblem(f,u₀=u₀))
end

"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=.5u-v``"
function poissonProblemExample_birthdeathinteractingsystem()
  f₁(u,x)  = ones(size(x,1)) - .5u[:,1]
  f₂(u,x)  = .5u[:,1] -   u[:,2]
  f(u,x) = [f₁(u,x) f₂(u,x)]
  u₀(x) = ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
  return(PoissonProblem(f,u₀=u₀))
end

## Stokes Examples

"Example problem for a homogeneous stationary Stokes equation."
function homogeneousStokesExample(;C=0)
  f₁(x,y)   = zeros(x)
  f₂(x,y)   = zeros(x)
  uanalytic(x,y) = 20x.*y.^3
  vanalytic(x,y) = 5x.^4 - 5y.^4
  panalytic(x,y) = 60x.^2 .*y - 20y.^3 + C
  g(x,y)    = 0
  return(StokesProblem(f₁,f₂,g,uanalytic,vanalytic,panalytic))
end

"Example problme for solving the trivial stationary Stokes equation."
function dirichletzeroStokesExample(;C=0)
  f₁(x,y)   = zeros(x)
  f₂(x,y)   = zeros(x)
  uanalytic(x,y) = 0
  vanalytic(x,y) = 0
  panalytic(x,y) = 0
  g(x,y)    = 0
  return(StokesProblem(f₁,f₂,g,uanalytic,vanalytic,panalytic))
end
