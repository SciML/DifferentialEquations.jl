##Boundary Setting Tests

using DifferentialEquations

mesh = meshExample_bunny()

findboundary(mesh,[1 4 5])
setboundary(mesh,"Robin")
