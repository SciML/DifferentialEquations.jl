srand(100)

### ODE Examples

# Linear ODE
f = (t,u) -> (1.01*u)
analytic = (t,u₀) -> u₀*exp(1.01*t)
"""Linear ODE on Float64"""
prob_ode_linear = ODEProblem(f,1/2,analytic=analytic)

const linear_bigα = parse(BigFloat,"1.01")
f = (t,u) -> (linear_bigα*u)
analytic = (t,u₀) -> u₀*exp(linear_bigα*t)
"""Linear ODE on Float64"""
prob_ode_bigfloatlinear = ODEProblem(f,parse(BigFloat,"0.5"),analytic=analytic)

f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = 1.01*u[i]
  end
end
analytic = (t,u₀) -> u₀*exp(1.01*t)
"""2D Linear ODE, standard 4x2"""
prob_ode_2Dlinear = ODEProblem(f,rand(4,2),analytic=analytic)
"""2D Linear ODE, 100x100"""
prob_ode_large2Dlinear = ODEProblem(f,rand(100,100),analytic=analytic)
"""2D Linaer ODE, bigfloats"""
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
prob_ode_bigfloat2Dlinear = ODEProblem(f,map(BigFloat,rand(4,2)).*ones(4,2)/2,analytic=analytic)
f = (t,u) -> 1.01*u
"""2D Linear ODE, not in place."""
prob_ode_2Dlinear_notinplace = ODEProblem(f,rand(4,2),analytic=analytic)

#Van der Pol Equations

f = (t,u,du) -> begin
  du[1] = ((1-u[2].^2)*u[1] - u[2])
  du[2] = u[1]
end
"""Van der Pol Equations. For non-stiff version"""
prob_ode_vanderpol = ODEProblem(f,[0;sqrt(3)])
f = (t,u,du) -> begin
  du[1] = ((1-u[2].^2)*u[1] - u[2])*1e6
  du[2] = u[1]
end
"""Van der Pol Equations. For difficult version."""
prob_ode_vanderpol_stiff = ODEProblem(f,[0;sqrt(3)])

# Lorenz Attractor

f = (t,u,du) -> begin
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - 8/3u[3]
end
"""Lorenz Attractor"""
prob_ode_lorenz = ODEProblem(f,ones(3))

# ROBER

f = (t,u,du) -> begin
  du[1] = -0.04*u[1]+1e4*u[2]*u[3]
  du[2] =  0.04*u[1]-3e7*(u[2])^2-1e4*u[2]*u[3]
  du[3] =  3e7*(u[2])^2
  nothing
end
"""
The Robertson biochemical reactions

http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Single/DemoRobertson/demorobertson.pdf

Or from Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 129

Usually solved on [0,1e11]
"""
prob_ode_rober = ODEProblem(f,[1.0;0.0;0.0])

# Three Body
const threebody_μ = parse(BigFloat,"0.012277471"); const threebody_μ′ = 1 - threebody_μ

f = (t,u,du) -> begin
  # 1 = y₁
  # 2 = y₂
  # 3 = y₁'
  # 4 = y₂'
  D₁ = ((u[1]+threebody_μ)^2 + u[2]^2)^(3/2)
  D₂ = ((u[1]-threebody_μ′)^2 + u[2]^2)^(3/2)
  du[1] = u[3]
  du[2] = u[4]
  du[3] = u[1] + 2u[4] - threebody_μ′*(u[1]+threebody_μ)/D₁ - threebody_μ*(u[1]-threebody_μ′)/D₂
  du[4] = u[2] - 2u[3] - threebody_μ′*u[2]/D₁ - threebody_μ*u[2]/D₂
end
"""
From Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 129

Usually solved on t₀ = 0.0; T = parse(BigFloat,"17.0652165601579625588917206249")
Periodic with that setup
"""
prob_ode_threebody = ODEProblem(f,[0.994, 0.0, 0.0, parse(BigFloat,"-2.00158510637908252240537862224")])

# Rigid Body Equations

f = (t,u,du) -> begin
  du[1]  = -2u[2]*u[3]
  du[2]  = 1.25*u[1]*u[3]
  du[3]  = -.5u[1]*u[2]
end

"""
Rigid Body Equations

From Solving Differential Equations in R by Karline Soetaert

or Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 244

Usually solved from 0 to 20. Periodic at 10 and 20.
"""
prob_ode_rigidbody = ODEProblem(f,[1.0,0.0,0.9])

# Pleiades Problem

f = (t,u,du) -> begin
  x = view(u,1:7)   # x
  y = view(u,8:14)  # y
  v = view(u,15:21) # x′
  w = view(u,22:28) # y′
  du[1:7] .= v
  du[8:14].= w
  du[14:21]=zero(u)
  for i=1:7,j=1:7
    if i != j
      r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
      du[14+i] += j*(x[j] - x[i])/r
      du[21+i] += j*(y[j] - y[i])/r
    end
  end
end
"""
Pleides Problem

From Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 244

Usually solved from 0 to 3.
"""
prob_ode_pleides = ODEProblem(f,[3.0,3.0,-1.0,-3.0,2.0,-2.0,2.0,3.0,-3.0,2.0,0,0,-4.0,4.0,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])

### SDE Examples

f = (t,u) -> 1.01*u
σ = (t,u) -> 0.87*u
analytic = (t,u₀,W) -> u₀.*exp(0.63155*t+0.87*W)
"""Example problem with solution ``u(t,W)=u₀*exp((α-(β^2)/2)*t+β*W)``"""
prob_sde_linear = SDEProblem(f,σ,1/2,analytic=analytic)

f = (t,u,du) -> begin
  for i = 1:length(u)
    du[i] = 1.01*u[i]
  end
end
σ = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = .87*u[i]
  end
end
"""Example problem of 8 linear SDEs (as a 4x2 matrix) with solution ``u(t,W)=u₀*exp((α-(β^2)/2)*t+β*W)``"""
prob_sde_2Dlinear = SDEProblem(f,σ,ones(4,2)/2,analytic=analytic)


f = (t,u) -> -.25*u*(1-u^2)
σ = (t,u) -> .5*(1-u^2)
analytic = (t,u₀,W) -> ((1+u₀).*exp(W)+u₀-1)./((1+u₀).*exp(W)+1-u₀)
"""Example problem with solution ``u(t,W)=((1+u₀)*exp(W)+u₀-1)./((1+u₀)*exp(W)+1-u₀)``"""
prob_sde_cubic = SDEProblem(f,σ,1/2,analytic=analytic)

f = (t,u) -> -0.01*sin(u).*cos(u).^3
σ = (t,u) -> 0.1*cos(u).^2
analytic = (t,u₀,W) -> atan(0.1*W + tan(u₀))
"""Example problem with solution ``u(t,W)=atan(0.1*W + tan(u₀))``"""
prob_sde_wave = SDEProblem(f,σ,1.,analytic=analytic)

const sde_wave_α = 0.1
const sde_wave_β = 0.05
f = (t,u) -> sde_wave_β./sqrt(1+t) - u./(2*(1+t))
σ = (t,u) -> sde_wave_α*sde_wave_β./sqrt(1+t)
analytic = (t,u₀,W) -> u₀./sqrt(1+t) + sde_wave_β*(t+sde_wave_α*W)./sqrt(1+t)
"""Example additive noise problem with solution ``u₀./sqrt(1+t) + β*(t+α*W)./sqrt(1+t)``"""
prob_sde_additive = SDEProblem(f,σ,1.,analytic=analytic)

const sde_wave_αvec = [0.1;0.1;0.1;0.1]
const sde_wave_βvec = [0.5;0.25;0.125;0.1115]
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = sde_wave_βvec[i]/sqrt(1+t) - u[i]/(2*(1+t))
  end
end

σ = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = sde_wave_αvec[i]*sde_wave_βvec[i]/sqrt(1+t)
  end
end
analytic = (t,u₀,W) -> u₀./sqrt(1+t) + sde_wave_βvec.*(t+sde_wave_αvec.*W)./sqrt(1+t)
"""Multiple Ito dimension extension of additiveSDEExample"""
prob_sde_additivesystem = SDEProblem(f,σ,[1.;1.;1.;1.],analytic=analytic)


f = (t,u,du) -> begin
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - 8/3u[3]
end
σ = (t,u) -> 3.0 #Additive
"""Lorenz Attractor with additive noise"""
prob_sde_lorenz = SDEProblem(f,σ,ones(3))

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
  function f(t,y)
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

  function σ1(t,y)
    dσ = zeros(19)
    dσ[1] = noiseLevel*1.5y[1]
    dσ[18]= noiseLevel*6y[18]
    return(dσ)
  end

  function σ2(t,y)
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

analytic_moving(t,x) = 0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x[:,1]-t+0.5).^2 + (x[:,2]-t+0.5).^2))
Du = (t,x) -> -50[analytic_moving(t,x).*(0.5-t+x[:,1])  analytic_moving(t,x).*(0.5-t+x[:,2])]
f = (t,x) -> (-5).*exp((-25).*((3/2)+6.*t.^2+x[:,1]+x[:,1].^2+x[:,2]+x[:,2].^2+(-2).*t.*(3+x[:,1]+
  x[:,2]))).*((-20)+(-100).*t.^2+(-49).*x[:,1]+(-50).*x[:,1].^2+(-49).*x[:,2]+(-50).*
  x[:,2].^2+2.*t.*(47+50.*x[:,1]+50.*x[:,2])+exp(25.*(1+(-2).*t).^2).*(22+
  100.*t.^2+49.*x[:,1]+50.*x[:,1].^2+49.*x[:,2]+50.*x[:,2].^2+(-2).*t.*(49+50.*x[:,1]+50.*x[:,2])))
"Example problem with solution: ``u(x,y,t)=0.1*(1-exp(-100*(t-0.5).^2)).*exp(-25((x-t+0.5).^2 + (y-t+0.5).^2))``"
prob_femheat_moving = HeatProblem(analytic_moving,Du,f)



analytic_diffuse(t,x) = exp(-10((x[:,1]-.5).^2 + (x[:,2]-.5).^2 )-t)
f = (t,x) -> exp(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
Du = (t,x) -> -20[analytic_diffuse(t,x).*(x[:,1]-.5) analytic_diffuse(t,x).*(x[:,2]-.5)]
"Example problem with solution: ``u(x,y,t)=exp(-10((x-.5).^2 + (y-.5).^2 )-t)``"
prob_femheat_diffuse = HeatProblem(analytic_diffuse,Du,f)


f = (t,x)  -> zeros(size(x,1))
u₀ = (x) -> float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6)) #Only mass at middle of (0,1)^2
"Example problem which starts with 1 at (0.5,0.5) and solves with ``f=gD=0``"
prob_femheat_pure = HeatProblem(u₀,f)


f = (t,x,u) -> ones(size(x,1)) - .5u
u₀ = (x) -> zeros(size(x,1))
"Example problem which starts with 0 and solves with ``f(u)=1-u/2``"
prob_femheat_birthdeath = HeatProblem(u₀,f)


f = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]   ones(size(x,1))-u[:,2]]
u₀ = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=1-v``"
prob_femheat_birthdeathsystem = HeatProblem(u₀,f)

f = (t,x,u)  -> [zeros(size(x,1))    zeros(size(x,1))]
u₀ = (x) -> [float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6)) float((abs(x[:,1]-.5) .< 1e-6) & (abs(x[:,2]-.5) .< 1e-6))]  # size (x,2), 2 meaning 2 variables
"Example problem which solves the homogeneous Heat equation with all mass starting at (1/2,1/2) with two different diffusion constants."
prob_femheat_diffusionconstants = HeatProblem(u₀,f,D=[.01 .001])

f  = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]     .5u[:,1]-u[:,2]]
u₀ = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables
"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=.5u-v``"
prob_femheat_birthdeathinteractingsystem = HeatProblem(u₀,f)

"Example problem which solves the Gray-Scott equations with quasi-random initial conditions"
function heatProblemExample_grayscott(;ρ=.03,k=.062,D=[1e-3 .5e-3])
  f₁(t,x,u)  = + u[:,1].*u[:,2].*u[:,2] + ρ*(1-u[:,2])
  f₂(t,x,u)  = u[:,1].*u[:,2].*u[:,2] -(ρ+k).*u[:,2]
  f(t,x,u) = [f₁(t,x,u) f₂(t,x,u)]
  u₀(x) = [ones(size(x,1))+rand(size(x,1)) .25.*float(((.2.<x[:,1].<.6) &
          (.2.<x[:,2].<.6)) | ((.85.<x[:,1]) & (.85.<x[:,2])))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

"Example problem which solves the Gierer-Meinhardt equations wtih quasi-random initial perturbations."
function heatProblemExample_gierermeinhardt(;a=1,α=1,D=[0.01 1.0],ubar=1,vbar=0,β=10,startNoise=0.01)
  f₁(t,x,u)  = a*u[:,1].*u[:,1]./u[:,2] + ubar - α*u[:,1]
  f₂(t,x,u)  = a*u[:,1].*u[:,1] + vbar -β.*u[:,2]
  f(t,x,u) = [f₁(t,x,u) f₂(t,x,u)]
  uss = (ubar +β)/α
  vss = (α/β)*uss.^2
  u₀(x) = [uss*ones(size(x,1))+startNoise*rand(size(x,1)) vss*ones(size(x,1))] # size (x,2), 2 meaning 2 variables
  return(HeatProblem(u₀,f,D=D))
end

f = (t,x,u)  -> ones(size(x,1)) - .5u
u₀ = (x) -> zeros(size(x,1))
σ = (t,x,u) -> 1u.^2
"Example problem which starts with 0 and solves with ``f(u)=1-u/2`` with noise ``σ(u)=10u^2``"
prob_femheat_stochasticbirthdeath = HeatProblem(u₀,f,σ=σ)

## Poisson

f = (x) -> sin(2π.*x[:,1]).*cos(2π.*x[:,2])
analytic = (x) -> sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
Du = (x) -> [cos(2*pi.*x[:,1]).*cos(2*pi.*x[:,2])./(4*pi) -sin(2π.*x[:,1]).*sin(2π.*x[:,2])./(4π)]
"Example problem with solution: ``u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)``"
prob_poisson_wave = PoissonProblem(f,analytic,Du)

σ = (x) -> 5 #Additive noise
"Example problem with deterministic solution: ``u(x,y)= sin(2π.*x).*cos(2π.*y)/(8π*π)``"
prob_poisson_noisywave = PoissonProblem(f,analytic,Du,σ=σ)

f = (x,u) -> ones(size(x,1)) - .5u
"Example problem for nonlinear Poisson equation. Uses ``f(u)=1-u/2``."
prob_poisson_birthdeath = PoissonProblem(f,numvars=1)

f  = (x,u) -> [ones(size(x,1))-.5u[:,1]     ones(size(x,1))-u[:,2]]
u₀ = (x) -> .5*ones(size(x,1),2) # size (x,2), 2 meaning 2 variables
"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=1-v``"
prob_poisson_birthdeathsystem = PoissonProblem(f,u₀=u₀)

f  = (x,u) -> [ones(size(x,1))-.5u[:,1]     .5u[:,1]-u[:,2]]
u₀ = (x) -> ones(size(x,1),2).*[.5 .5] # size (x,2), 2 meaning 2 variables

"Example problem which starts with 1/2 and solves the system ``f(u)=1-u/2`` and ``f(v)=.5u-v``"
prob_poisson_birthdeathinteractingsystem = PoissonProblem(f,u₀=u₀)

## Stokes Examples

f₁ = (x,y)   -> zeros(x)
f₂ = (x,y)   -> zeros(x)
uanalytic = (x,y) -> 20x.*y.^3
vanalytic = (x,y) -> 5x.^4 - 5y.^4
panalytic = (x,y) -> 60x.^2 .*y - 20y.^3# + C
g = (x,y)    -> 0

"Example problem for a homogeneous stationary Stokes equation."
prob_stokes_homogenous = StokesProblem(f₁,f₂,g,uanalytic,vanalytic,panalytic)

f₁ = (x,y)   -> zeros(x)
f₂ = (x,y)   -> zeros(x)
uanalytic = (x,y) -> 0
vanalytic = (x,y) -> 0
panalytic = (x,y) -> 0
g = (x,y)    -> 0
"Example problme for solving the trivial stationary Stokes equation."
prob_stokes_dirichletzero = StokesProblem(f₁,f₂,g,uanalytic,vanalytic,panalytic)
