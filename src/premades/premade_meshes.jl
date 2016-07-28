
###Example Meshes
pkgdir = Pkg.dir("DifferentialEquations")
meshes_location = "src/premades/premade_meshes.jld"
"meshExample_bunny() : Returns a 3D SimpleMesh."
meshExample_bunny() = Main.JLD.load("$pkgdir/$meshes_location","bunny")

"meshExample_flowpastcylindermesh() : Returns a 2D SimpleMesh."
meshExample_flowpastcylindermesh() = Main.JLD.load("$pkgdir/$meshes_location","flowpastcylindermesh")

"meshExample_lakemesh() : Returns a 2D SimpleMesh."
meshExample_lakemesh() = Main.JLD.load("$pkgdir/$meshes_location","lakemesh")

"meshExample_Lshapemesh() : Returns a 2D SimpleMesh."
meshExample_Lshapemesh() = Main.JLD.load("$pkgdir/$meshes_location","Lshapemesh")

"meshExample_Lshapeunstructure() : Returns a 2D SimpleMesh."
meshExample_Lshapeunstructure() = Main.JLD.load("$pkgdir/$meshes_location","Lshapeunstructure")

"meshExample_oilpump() : Returns a 3D SimpleMesh."
meshExample_oilpump() = Main.JLD.load("$pkgdir/$meshes_location","oilpump")

"meshExample_wavymesh() : Returns a 2D SimpleMesh."
meshExample_wavymesh() = Main.JLD.load("$pkgdir/$meshes_location","wavymesh")

"meshExample_wavyperturbmesh() : Returns a 3D SimpleMesh."
meshExample_wavyperturbmesh() = Main.JLD.load("$pkgdir/$meshes_location","wavyperturbmesh")
