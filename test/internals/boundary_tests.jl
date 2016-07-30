##Boundary Setting Tests

using DifferentialEquations, JLD

mesh = meshExample_lakemesh()

findboundary(mesh,ones(size(vec(mesh.elem))))
setboundary(mesh,:robin)

true
