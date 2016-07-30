"""
GSu!(u,f₁,Δxs,p,ugD,grids,ux,uy)

Performs a Gauss-Seidel iteration on u.
"""
function GSu!(u,f₁,Δxs,p,ugD,grids,ux,uy)
  inti = 2:size(u,1)-1; intj = 2:size(u,2)-1 # Coordinates for interior
  fh₁ = Δxs[1]*Δxs[2]*f₁(ux,uy) # Get f values plus ghosts, include left/right boundary

  # Now solve, dirichlet boundary conditions pre-imposed on verticals
  for j in intj #Top boundary
    bval = fh₁[1,j] + (8/3)*ugD(grids[1][1,j],grids[2][1,j])
    u[1,j]  = (bval + (4/3)*u[2,j] + u[1,j-1] + u[1,j+1] - Δxs[1]*(p[1,j]-p[1,j-1]))/6
  end
  # Solve on interior
  for j in intj, i in inti #Interior,
    u[i,j] = (fh₁[i,j] + u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - Δxs[1]*(p[i,j]-p[i,j-1]))/4
  end
  for j in intj #Bottom boundary
    bval = fh₁[end,j] + (8/3)*ugD(grids[1][end,j],grids[2][end,j])
    u[end,j]= (bval + (4/3)*u[end-1,j] + u[end,j-1] + u[end,j+1] -Δxs[1]*(p[end,j]-p[end,j-1]))/6
  end
end

"""
GSv!(v,f₂,Δxs,p,vgD,grids,vx,vy)

Performs a Gauss-Seidel iteration on v.
"""
function GSv!(v,f₂,Δxs,p,vgD,grids,vx,vy)
  inti = 2:size(v,1)-1; intj = 2:size(v,2)-1 # Coordinates for interior
  fh₂   = Δxs[1]*Δxs[2]*f₂(vx,vy) # Get f values plus ghosts, include top/bottom boundary

  # Now solve, dirichlet boundary conditions pre-imposed on horizontals
  for i in inti #Left boundary
    bval = fh₂[i,1] + (8/3)*vgD(grids[1][i,1],grids[2][i,1])
    v[i,1]  = (bval + (4/3)*v[i,2] + v[i-1,1] + v[i+1,1] - Δxs[2]*(p[i,1]-p[i-1,1]))/6
  end
  # Solve on interior
  for j in intj,i in inti  #Interior
    v[i,j] = (fh₂[i,j] + v[i-1,j] + v[i+1,j] + v[i,j-1] + v[i,j+1] - Δxs[2]*(p[i,j]-p[i-1,j]))/4
  end
  for i in inti #Right boundary
    bval = fh₂[i,end] + (8/3)*vgD(grids[1][i,end],grids[2][i,end])
    v[i,end]  = (bval + (4/3)*v[i,end-1] + v[i-1,end] + v[i+1,end] - Δxs[2]*(p[i,end]-p[i-1,end]))/6
  end
end

"""
calc_rp!(rp,u,v,Δxs,g,px,py)

Calculates the rp from the u and v's.
"""
function calc_rp!(rp,u,v,Δxs,g,px,py)
  for j=1:size(rp,2),i=1:size(rp,1)
    rp[i,j] = Δxs[1]*Δxs[2]*g(px[i,j],py[i,j]) - (u[i,j+1]-u[i,j])*Δxs[2] - (v[i+1,j]-v[i,j])*Δxs[1]
  end
end

"""
uzawa_p!(p,u,v,Δxs,g,px,py)

Solves for p from u and v using an Uzawa update.
"""
function uzawa_p!(p,u,v,Δxs,g,px,py)
  for j=1:size(p,2),i=1:size(p,1)
    p[i,j] = p[i,j] + 0.01(-(u[i,j+1]-u[i,j])/Δxs[1] - (v[i+1,j]-v[i,j])/Δxs[2] + g(px[i,j],py[i,j]))
  end
end

"""
GSδq!(δq,rp,Δxs)

Performs a Gauss-Seidel iteration for δq.
"""
function GSδq!(δq,rp,Δxs)
  for j=2:size(δq,2)-1, i=2:size(δq,1)-1 # Interior
    δq[i,j] = (rp[i,j] + δq[i-1,j] + δq[i+1,j] + δq[i,j-1] + δq[i,j+1])/4
  end
  for j=2:size(δq,2)-1
    i=1
    δq[i,j] = (rp[i,j] + δq[i+1,j] + δq[i,j-1] + δq[i,j+1])/3
    i=size(δq,1)
    δq[i,j] = (rp[i,j] + δq[i-1,j]  + δq[i,j-1] + δq[i,j+1])/3
  end
  for i=2:size(δq,1)-1
    j=1
    δq[i,j] = (rp[i,j] + δq[i-1,j] + δq[i+1,j] + δq[i,j+1])/3
    j=size(δq,2)
    δq[i,j] = (rp[i,j] + δq[i-1,j] + δq[i+1,j] + δq[i,j-1])/3
  end
  #Corners
  i=1; j=1
  δq[i,j] = (rp[i,j] + δq[i+1,j] + δq[i,j+1])/2
  i=size(δq,1)
  δq[i,j] = (rp[i,j] + δq[i-1,j] + δq[i,j+1])/2
  i=1; j=size(δq,2)
  δq[i,j] = (rp[i,j] + δq[i+1,j] + δq[i,j-1])/2
  i=size(δq,1)
  δq[i,j] = (rp[i,j] + δq[i-1,j] + δq[i,j-1])/2
end

"""
update_u!(u,δq,Δxs)

Updates u given δq
"""
function update_u!(u,δq,Δxs)
  for i=1:size(u,1),j=2:size(u,2)-1
    u[i,j] = u[i,j] - (δq[i,j]-δq[i,j-1])/Δxs[2]
  end
end

"""
update_v!(v,δq,Δxs)

Updates v given δq
"""
function update_v!(v,δq,Δxs)
  for i=2:size(v,1)-1,j=1:size(v,2)
    v[i,j] = v[i,j] - (δq[i,j]-δq[i-1,j])/Δxs[1]
  end
end

"""
update_p!(p,δq,Δxs)

Updates p given δq
"""
function update_p!(p,δq,Δxs)
  for j=2:size(p,2)-1, i=2:size(p,1)-1 #Interior
    p[i,j] = p[i,j] - (4δq[i,j]-δq[i-1,j]-δq[i+1,j]-δq[i,j-1]-δq[i,j+1])/(Δxs[1]*Δxs[2])
  end
  for j=2:size(p,2)-1
    i=1
    p[i,j] = p[i,j] - (3δq[i,j]-δq[i+1,j]-δq[i,j-1]-δq[i,j+1])/(Δxs[1]*Δxs[2])
    i=size(p,1)
    p[i,j] = p[i,j] - (3δq[i,j]-δq[i-1,j]-δq[i,j-1]-δq[i,j+1])/(Δxs[1]*Δxs[2])
  end
  for i=2:size(p,1)-1
    j=1
    p[i,j] = p[i,j] - (3δq[i,j]-δq[i-1,j]-δq[i+1,j]-δq[i,j+1])/(Δxs[1]*Δxs[2])
    j=size(p,2)
    p[i,j] = p[i,j] - (3δq[i,j]-δq[i-1,j]-δq[i+1,j]-δq[i,j-1])/(Δxs[1]*Δxs[2])
  end
  #Corners
  i=1; j=1
  p[i,j] = p[i,j] - (2δq[i,j]-δq[i+1,j]-δq[i,j+1])/(Δxs[1]*Δxs[2])
  i=size(p,1)
  p[i,j] = p[i,j] - (2δq[i,j]-δq[i-1,j]-δq[i,j+1])/(Δxs[1]*Δxs[2])
  i=1;j=size(p,2)
  p[i,j] = p[i,j] - (2δq[i,j]-δq[i+1,j]-δq[i,j-1])/(Δxs[1]*Δxs[2])
  i=size(p,1)
  p[i,j] = p[i,j] - (2δq[i,j]-δq[i-1,j]-δq[i,j-1])/(Δxs[1]*Δxs[2])
end

"""
stokes_restriction(u,v,p,Δxs,grids,mins,maxs,ugD,vgD)

Restricts the Stokes problem to the coarsegrid.
"""
function stokes_restriction(u,v,p,Δxs,grids,mins,maxs,ugD,vgD)
  n = size(p,1); m = size(p,2)
  newu = Array{Float64}(n÷2,m÷2+1)
  newv = Array{Float64}(n÷2+1,m÷2)
  newp = Array{Float64}(n÷2,m÷2)
  if n%2!==0 || m%2!==0
    error("Warning, size must be divisible by two. Aborting")
  end
  #Calculate new grids
  Δxs = Δxs*2
  grids = meshgrid(mins[1]:Δxs[1]:maxs[1],mins[2]:Δxs[2]:maxs[2])
  ux = grids[1][1:end-1,:]
  uy = (grids[2]+Δxs[2]/2)[1:end-1,:]
  vx = (grids[1]+Δxs[1]/2)[:,1:end-1]
  vy = grids[2][:,1:end-1]
  px = grids[1][1:end-1,1:end-1]+Δxs[1]/2
  py = grids[2][1:end-1,1:end-1]+Δxs[2]/2
  #Restrict u interior
  for j=2:m÷2,i=1:n÷2
    newu[i,j] = (u[2i-1,2j-2]+2u[2i-1,2j-1]+u[2i-1,2j]+u[2i,2j-2]+2u[2i,2j-1]+u[2i,2j])/8
  end
  #Restrict v interior
  for j=1:m÷2,i=2:n÷2
    newv[i,j] = (v[2i-2,2j-1]+v[2i-2,2j]+2v[2i-1,2j-1]+2v[2i-1,2j]+v[2i,2j-1]+v[2i,2j])/8
  end
  #Set u on boundary
  newu[:,1]  =ugD(ux[:,1]  ,uy[:,1])
  newu[:,end]=ugD(ux[:,end],uy[:,end])
  newv[1,:]  =vgD(vx[1,:]  ,vy[1,:])
  newv[end,:]=vgD(vx[end,:],vy[end,:])
  #Restrict p
  for j = 1:m÷2, i = 1:n÷2
    newp[i,j] = (p[2i,2j]+p[2i-1,2j]+p[2i,2j-1]+p[2i-1,2j-1])/4
  end
  return(newu,newv,newp,Δxs,grids,ux,uy,vx,vy,px,py)
end

"""
stokes_prolongation(u,v,p,Δxs,grids,mins,maxs,ugD,vgD)

Prolongates the Stokes problem to the fine grid
"""
function stokes_prolongation(u,v,p,Δxs,grids,mins,maxs,ugD,vgD)
  n = size(p,1); m = size(p,2)
  newu = Array{Float64}(2n,2m+1)
  newv = Array{Float64}(2n+1,2m)
  newp = Array{Float64}(2n,2m)
  if n%2!==0 || m%2!==0
    error("Warning, size must be divisible by two. Aborting")
  end
  #Calculate new grids
  Δxs = Δxs/2
  grids = meshgrid(mins[1]:Δxs[1]:maxs[1],mins[2]:Δxs[2]:maxs[2])
  ux = grids[1][1:end-1,:]
  uy = (grids[2]+Δxs[2]/2)[1:end-1,:]
  vx = (grids[1]+Δxs[1]/2)[:,1:end-1]
  vy = grids[2][:,1:end-1]
  px = grids[1][1:end-1,1:end-1]+Δxs[1]/2
  py = grids[2][1:end-1,1:end-1]+Δxs[2]/2
  #Prolongate u interior
  for j=1:m-1,i=1:n-1
    u₁ = u[i,j]; u₂ = u[i+1,j]; u₃ = u[i,j+1]; u₄ = u[i+1,j+1]
    newu[2i,2j]    = (3/4)*u₁ + (1/4)*u₂ #u₅
    newu[2i+1,2j]  = (3/4)*u₂ + (1/4)*u₁ #u₆
    newu[2i,2j+2]  = (3/4)*u₃ + (1/4)*u₄ #u₉
    newu[2i+1,2j+2]= (3/4)*u₄ + (1/4)*u₃ #u₁₀
    newu[2i,2j+1]  = (1/2)*(newu[2i,2j]+newu[2i,2j+2]) #u₇
    newu[2i+1,2j+1]= (1/2)*(newu[2i+1,2j]+newu[2i+1,2j+2]) #u₈
  end
  #Prolongate v interior
  for j=1:m-1,i=1:n-1
    v₁ = v[i,j]; v₂ = v[i+1,j]; v₃ = v[i,j+1]; v₄ = v[i+1,j+1]
    newv[2i,2j]    = (3/4)*v₁ + (1/4)*v₃ #v₅
    newv[2i+2,2j]  = (3/4)*v₂ + (1/4)*v₄ #v₇
    newv[2i,2j+1]  = (3/4)*v₃ + (1/4)*v₁ #v₈
    newv[2i+2,2j+1]= (3/4)*v₄ + (1/4)*v₂ #v₁₀
    newv[2i+1,2j]  = (1/2)*(newv[2i,2j]+newv[2i+2,2j]) #v₆
    newv[2i+1,2j+1]= (1/2)*(newv[2i,2j+1]+newv[2i+2,2j+1]) #v₉
  end
  #Get Top/Bottom Boundaries
  for j=1:m-1
    newu[1,2j]    =(1/2)*(u[1,j]  +ugD(grids[1][1,2j],grids[2][1,2j]))
    newu[1,2j+2]  =(1/2)*(u[1,j+1]+ugD(grids[1][1,2j+2],grids[2][1,2j+2]))
    newu[1,2j+1]  =(1/4)*(u[1,j]+u[1,j+1]+2ugD(grids[1][1,2j+1],grids[2][1,2j+1]))
    newu[end,2j]  =(1/2)*(u[end,j]  +ugD(grids[1][end,2j],grids[2][end,2j]))
    newu[end,2j+2]=(1/2)*(u[end,j+1]+ugD(grids[1][end,2j+2],grids[2][end,2j+2]))
    newu[end,2j+1]=(1/4)*(u[end,j]+u[end,j+1]+2ugD(grids[1][end,2j+1],grids[2][end,2j+1]))
  end
  #Get Left/Right Boundaries
  for i=1:n-1
    newv[2i,1]    =(1/2)*(v[i,1]  +vgD(grids[1][2i,1],grids[2][2i,1]))
    newv[2i+2,1]  =(1/2)*(v[i+1,1]+vgD(grids[1][2i+2,1],grids[2][2i+2,1]))
    newv[2i+1,1]  =(1/4)*(v[i,1]+v[i+1,1]+2vgD(grids[1][2i+1,1],grids[2][2i+1,1]))
    newv[2i,end]  =(1/2)*(v[i,end]  +vgD(grids[1][2i,end],grids[2][2i,end]))
    newv[2i+2,end]=(1/2)*(v[i+1,end]+vgD(grids[1][2i+2,end],grids[2][2i+2,end]))
    newv[2i+1,end]=(1/4)*(v[i,end]+v[i+1,end]+2vgD(grids[1][2i+1,end],grids[2][2i+1,end]))
  end
  #Set u on boundary

  newu[:,1]  =ugD(ux[:,1]  ,uy[:,1])
  newu[:,end]=ugD(ux[:,end],uy[:,end])
  newv[1,:]  =vgD(vx[1,:]  ,vy[1,:])
  newv[end,:]=vgD(vx[end,:],vy[end,:])

  #Prolongate p
  for j = 1:2m, i = 1:2n
    newp[i,j] = p[(i-1)÷2+1,(j-1)÷2+1]
  end
  return(newu,newv,newp,Δxs,grids,ux,uy,vx,vy,px,py)
end

@def dgs begin
  for j = 1:gsiters
    GSu!(u,f₁,Δxs,p,ugD,grids,ux,uy) # Inplace u => uhalf
    GSv!(v,f₂,Δxs,p,vgD,grids,vx,vy) # Inplace v => vhalf
  end
  calc_rp!(rp,u,v,Δxs,g,px,py)
  for j = 1:gsiters
    GSδq!(δq,rp,Δxs)
  end
  δq = -δq # Find out how to fix this. Too hacky
  update_u!(u,δq,Δxs)
  update_v!(v,δq,Δxs)
  update_p!(p,δq,Δxs)
end

"""
solve(prob::StokesProblem,mesh::FDMMesh)

Solves the given stationary Stokes problem on the given finite difference mesh.

### Keyword Arguments

* `converrors`: Whether to calculate all of the errors along the convergence. Default is true.
* `maxiters`: Maximum number of iterations before haulting. Default is 100.
* `alg`: The solver algorithm. Default is "dgs". Other option is "multigrid".
* `level`: The number of levels in the Multigrid. Default is 2.
* `smoothSteps`: The number of Gauss-Seidel iterations to do at each smoothing step. Default is 10.
* `coarseSteps`: The number of Gauss-Seidel iterations to do at the coarsegrid. Default is 40.
* `gsiters`: The number of Gauss-Seidel iterations to do at each step. Default is 20.
"""
function solve(prob::StokesProblem,mesh::FDMMesh;converrors=true,maxiters=100,alg=:DGS,levels=2,smoothSteps=10,coarseSteps=40,gsiters=20)
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

  if trueKnown
    uTrue = float(uanalytic(ux,uy))
    vTrue = float(vanalytic(vx,vy))
    pTrue = float(panalytic(px,py))
  else
    uTrue = nothing
    vTrue = nothing
    pTrue = nothing
  end

  # Impose boundary conditions, cut out ghost points
  u[:,1]   = ugD(ux[:,1]  ,uy[:,1])
  u[:,end] = ugD(ux[:,end],uy[:,end])
  v[1,:]   = vgD(vx[1,:]  ,vy[1,:])
  v[end,:] = vgD(vx[end,:],vy[end,:])

  #Note if Atom is loaded for progress
  atomloaded = isdefined(Main,:Atom)

  if converrors
    if !trueKnown
      error("True not known. Cannot calculate errors. Abandoning convergence calculations")
      converrors = false
    end
    uold = float(uanalytic(ux,uy))
    vold = float(vanalytic(vx,vy))
    pold = float(panalytic(px,py))
    converror_maxu    = Vector{Float64}(0)
    converror_maxv    = Vector{Float64}(0)
    converror_maxp    = Vector{Float64}(0)
    converror_l2u     = Vector{Float64}(0)
    converror_l2v     = Vector{Float64}(0)
    converror_l2p     = Vector{Float64}(0)
    converror_relmaxu = Vector{Float64}(0)
    converror_relmaxv = Vector{Float64}(0)
    converror_relmaxp = Vector{Float64}(0)
    converror_rell2u  = Vector{Float64}(0)
    converror_rell2v  = Vector{Float64}(0)
    converror_rell2p  = Vector{Float64}(0)
    converror_relresl2u  = Vector{Float64}(0)
    converror_relresl2v  = Vector{Float64}(0)
    converror_relresl2p  = Vector{Float64}(0)
  end
  for i = 1:maxiters
    if converrors
      uold[:] = u
      vold[:] = v
      pold[:] = p
    end
    if alg==:DGS
      @dgs
    elseif alg==:Multigrid
      j=1
      while j<levels
        j+=1
        for k=1:smoothSteps
          @dgs # Smooth
        end
        u,v,p,Δxs,grids,ux,uy,vx,vy,px,py = stokes_restriction(u,v,p,Δxs,grids,mins,maxs,ugD,vgD) # Restrict
        rp = zeros(p)
        δq = zeros(p)
      end
      for k=1:coarseSteps
        @dgs #Solve quasi-exact
      end
      j=1
      while j<levels
        j+=1
        u,v,p,Δxs,grids,ux,uy,vx,vy,px,py = stokes_prolongation(u,v,p,Δxs,grids,mins,maxs,ugD,vgD) #Prolongate
        rp = zeros(p)
        δq = zeros(p)
        for k=1:smoothSteps
          @dgs # Smooth
        end
      end
    end
    if converrors
      push!(converror_maxu,maximum(abs(u-uTrue)))
      push!(converror_maxv,maximum(abs(v-vTrue)))
      push!(converror_maxp,maximum(abs(p-pTrue)))
      push!(converror_l2u,norm(u-uTrue,2))
      push!(converror_l2v,norm(v-vTrue,2))
      push!(converror_l2p,norm(p-pTrue,2))
      push!(converror_relmaxu,maximum(abs(u-uTrue))/maximum(abs(u)))
      push!(converror_relmaxv,maximum(abs(v-vTrue))/maximum(abs(v)))
      push!(converror_relmaxp,maximum(abs(p-pTrue))/maximum(abs(p)))
      push!(converror_rell2u,norm(u-uTrue,2)/norm(u,2))
      push!(converror_rell2v,norm(v-vTrue,2)/norm(v,2))
      push!(converror_rell2p,norm(p-pTrue,2)/norm(p,2))
      push!(converror_relresl2u,norm(u-uold,2)/norm(u,2))
      push!(converror_relresl2v,norm(v-vold,2)/norm(v,2))
      push!(converror_relresl2p,norm(p-pold,2)/norm(p,2))
    end
    atomloaded ? Main.Atom.progress(i/maxiters) : nothing
  end

  #Generate and return solution type
  if trueKnown
    errors = Dict(:ul∞=>maximum(abs(u-uTrue)),:vl∞=>maximum(abs(v-vTrue)),
                  :pl∞=>maximum(abs(p-pTrue)),:ul2=>norm(u-uTrue,2),
                  :vl2=>norm(v-vTrue,2),:pl2=>norm(p-pTrue,2),
                  :rul∞=>maximum(abs(u-uTrue))/maximum(abs(u)),
                  :rvl∞=>maximum(abs(v-vTrue))/maximum(abs(v)),
                  :rpl∞=>maximum(abs(p-pTrue))/maximum(abs(p)),
                  :rul2=>norm(u-uTrue,2)/norm(u,2),
                  :rvl2=>norm(v-vTrue,2)/norm(v,2),
                  :rpl2=>norm(p-pTrue,2)/norm(p,2),
    )
  else
    error = nothing
  end
  if converrors
    converrors = Dict(:ul∞=>converror_maxu,:vl∞=>converror_maxv,
    :pl∞=>converror_maxp,:ul2=>converror_l2u,
    :vl2=>converror_l2v,:pl2=>converror_l2p,
    :rul∞=>converror_relmaxu,:rvl∞=>converror_relmaxv,
    :rpl∞=>converror_relmaxp,:rul2=>converror_rell2u,
    :rvl2=>converror_rell2v,:rpl2=>converror_rell2p,
    :rresul2=>converror_relresl2u,
    :rresvl2=>converror_relresl2v,:rrespl2=>converror_relresl2p,
    )
  else
    converrors = nothing
  end

  return(StokesSolution(u,v,p,uTrue,vTrue,pTrue,mesh,trueKnown;errors=errors,converrors=converrors))
end
