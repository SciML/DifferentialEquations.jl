##Boundary Setting Tests

using DifferentialEquations

mesh = meshExample_bunny()

findboundary(mesh,ones(size(vec(mesh.elem))))
setboundary(mesh,"Robin")

true
