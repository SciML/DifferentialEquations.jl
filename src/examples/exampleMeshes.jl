
###Example Meshes
pkgdir = Pkg.dir("DifferentialEquations")
"meshExample_bunny() : Returns a 3D SimpleMesh."
meshExample_bunny() = load("$pkgdir/src/examples/prebuiltMeshes.jld","bunny")

"meshExample_flowpastcylindermesh() : Returns a 2D SimpleMesh."
meshExample_flowpastcylindermesh() = load("$pkgdir/src/examples/prebuiltMeshes.jld","flowpastcylindermesh")

"meshExample_lakemesh() : Returns a 2D SimpleMesh."
meshExample_lakemesh() = load("$pkgdir/src/examples/prebuiltMeshes.jld","lakemesh")

"meshExample_Lshapemesh() : Returns a 2D SimpleMesh."
meshExample_Lshapemesh() = load("$pkgdir/src/examples/prebuiltMeshes.jld","Lshapemesh")

"meshExample_Lshapeunstructure() : Returns a 2D SimpleMesh."
meshExample_Lshapeunstructure() = load("$pkgdir/src/examples/prebuiltMeshes.jld","Lshapeunstructure")

"meshExample_oilpump() : Returns a 3D SimpleMesh."
meshExample_oilpump() = load("$pkgdir/src/examples/prebuiltMeshes.jld","oilpump")

"meshExample_wavymesh() : Returns a 2D SimpleMesh."
meshExample_wavymesh() = load("$pkgdir/src/examples/prebuiltMeshes.jld","wavymesh")

"meshExample_wavyperturbmesh() : Returns a 3D SimpleMesh."
meshExample_wavyperturbmesh() = load("$pkgdir/src/examples/prebuiltMeshes.jld","wavyperturbmesh")
