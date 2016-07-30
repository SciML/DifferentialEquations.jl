###Example Meshes

using JLD

mesh = meshExample_bunny()

mesh = meshExample_flowpastcylindermesh()

mesh = meshExample_lakemesh()

mesh = meshExample_Lshapemesh()

mesh = meshExample_Lshapeunstructure()

mesh = meshExample_oilpump()

mesh = meshExample_wavymesh()

mesh = meshExample_wavyperturbmesh()

TEST_PLOT && plot(mesh)

true
