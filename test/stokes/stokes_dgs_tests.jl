using DifferentialEquations, Parameters, Plots
Δx = 1//2^3 # Make this much smaller (1^2-5) for your own tests
mesh = FDMMesh(Δx,mins=[-1;-1],maxs=[1;1])
prob = dirichletzeroStokesExample()
prob = homogeneousStokesExample()

@unpack mesh: Δxs,grids,dims,gridSize,square,mins,maxs
u = zeros(gridSize[1]-1,gridSize[2])
v = zeros(gridSize[1],gridSize[2]-1)
p = zeros(gridSize[1]-1,gridSize[2]-1)
rp = zeros(p)
δq = zeros(p)

#Calculate grids for u,v, and p
ux = grids[1][1:end-1,:]
uy = (grids[2]+Δxs[2]/2)[1:end-1,:]
vx = (grids[1]+Δxs[1]/2)[:,1:end-1]
vy = grids[2][:,1:end-1]
px = grids[1][1:end-1,1:end-1]+Δxs[1]/2
py = grids[2][1:end-1,1:end-1]+Δxs[2]/2
#prob = dirichletzeroStokesExample()
@unpack prob: f₁,f₂,ugD,vgD,uanalytic,vanalytic,panalytic,g,trueKnown

u_analytic = float(uanalytic(ux,uy))
vTrue = float(vanalytic(vx,vy))
pTrue = float(panalytic(px,py))

# Impose boundary conditions, cut out ghost points
u[:,1]   = ugD(ux[:,1]  ,uy[:,1])
u[:,end] = ugD(ux[:,end],uy[:,end])
v[1,:]   = vgD(vx[1,:]  ,vy[1,:])
v[end,:] = vgD(vx[end,:],vy[end,:])

#Gauss-Seidel of Velocity Test
err1 = Vector{Float64}(0)
err2 = Vector{Float64}(0)
DifferentialEquations.GSu!(u,f₁,Δxs,pTrue,ugD,grids,ux,uy) # Inplace u => uhalf
DifferentialEquations.GSv!(v,f₂,Δxs,pTrue,vgD,grids,vx,vy) # Inplace v => vhalf
for j = 1:100
  DifferentialEquations.GSu!(u,f₁,Δxs,pTrue,ugD,grids,ux,uy) # Inplace u => uhalf
  DifferentialEquations.GSv!(v,f₂,Δxs,pTrue,vgD,grids,vx,vy) # Inplace v => vhalf
  DifferentialEquations.push!(err1,maximum(u-u_analytic))
  DifferentialEquations.push!(err2,maximum(v-vTrue))
end

#Should Converge
TEST_PLOT && plot(1:100,[err1 err2],yscale=:log10)

#Now test rp
err1 = Vector{Float64}(0)
val  = Vector{Float64}(0)

DifferentialEquations.calc_rp!(rp,u_analytic,vTrue,Δxs,g,px,py)

#Should be close to zero with true values.
maximum(rp) < 1e-2
