using DifferentialEquations

## Order 1

a = 0
b = 1
l = 1

f1(x,n) = x.^(2*n-1)
f2(x) = sin(x)
err = Array{Float64}(20,2)
for n=1:20
  λ,ω = quadpts1(n)
  nQuad = size(λ,1)
  t1 = 0
  t2 = 0
  for p = 1:nQuad
      px = λ[p,1]*a+ λ[p,2]*b
      t1 = t1 + ω[p]*f1(px,n)
      t2 = t2 + ω[p]*f2(px)
  end
  t1 = t1*l
  t2 = t2*l
  err[n,1] = abs(t1 - 1/(2*n))
  err[n,2] = abs(t2 - (- cos(1) - 1))
end
println("Quadrature Errors")
println(err)

## Order 2

#Functions to test integration on
f1(x,y,n) = x.^n + y.^n #Error should be exact
f2(x,y) = sin(x+y) #Error should decay with n
node = [0 0; 1 0; 0 1] # Reference Triangle
elem = [1 2 3]
area = 0.5
err = Array{Float64}(10,2)
for n = 1:10 #Integrate with each n
    λ,ω = quadpts(n)
    nQuad = size(λ,1)
    t1 = 0
    t2 = 0
    for p = 1:nQuad
        pxy = λ[p,1]*node[elem[:,1],:] +
              λ[p,2]*node[elem[:,2],:] +
              λ[p,3]*node[elem[:,3],:]
        t1 = t1 + ω[p]*f1(pxy[1],pxy[2],n)
        t2 = t2 + ω[p]*f2(pxy[1],pxy[2])
    end
    t1 = t1*area
    t2 = t2*area
    err[n,1] = abs(t1 - 2/((n+1)*(n+2)))
    err[n,2] = abs(t2 - (sin(1) - cos(1)))
end
println("Quadrature Errors")
println(err)

minimum(err[1:end-1,1] .< 1e-14) && minimum(diff(err[1:end-1,2]) .<0)
