## Test of Assembly Materials

using DifferentialEquations

Δx = 1//2^(5)
femMesh = notime_squaremesh([0 1 0 1],Δx,"Dirichlet")


A,M,area = assemblematrix(femMesh.node,femMesh.elem,lumpflag=false,K=4)

abs(M[2,1] - 40.6901e-5) < .001 #Non-diagonal
