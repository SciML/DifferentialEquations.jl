#Meshgrid Tests

using DifferentialEquations

bool1 = size(meshgrid(0:0.1:1),1) == 2
bool2 = size(meshgrid(0:0.1:1)[1]) == (11,11)
bool3 = meshgrid(0:0.1:1) == meshgrid(0:0.1:1,0:0.1:1)
bool4 = size(meshgrid(0:0.1:1,0:0.1:1,0:0.1:1)[1]) == (11,11,11)

bool1 && bool2 && bool3 && bool4
