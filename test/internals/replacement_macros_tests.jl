using DifferentialEquations

f = (t,u,du) -> begin
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end

g = @ode_define begin
  dx = σ*(y-x)
  dy = x*(ρ-z) - y
  dz = x*y - β*z
end σ=>10. ρ=>28. β=>(8/3)

du = zeros(3)
du2= zeros(3)
f(1.0,[1.0,2.0,3.0],du)
g(1.0,[1.0,2.0,3.0],du2)

f = @fem_define((t,x),(),exp(-t-5*(1-2x+2x.^2 - 2y +2y.^2)).*(-161 + 400*(x - x.^2 + y - y.^2)))
g = (t,x) -> exp(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
x = rand(4,2)

h = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]   ones(size(x,1))-u[:,2]]

l = @fem_define((t,x,u),(u,v),begin
  [ones(length(u))-α*u ones(length(v))-v]
end,α=>0.5)

f(1.0,x) == g(1.0,x)
h(1.0,x,x) == l(1.0,x,x)
du == du2
