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

du == du2
