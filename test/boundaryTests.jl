##Boundary Setting Tests

using DifferentialEquations

mesh = meshExample_lakemesh()

findboundary(mesh,ones(size(vec(mesh.elem))))
setboundary(mesh,"Robin")

true
