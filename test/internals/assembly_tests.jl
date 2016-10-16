## Test of Assembly Materials

using DifferentialEquations

Δx = 1//2^(5)
fem_mesh = notime_squaremesh([0 1 0 1],Δx,:dirichlet)


A,M,area = assemblematrix(fem_mesh.node,fem_mesh.elem,lumpflag=false,K=4)

abs.(M[2,1] - 40.6901e-5) < .001 #Non-diagonal
