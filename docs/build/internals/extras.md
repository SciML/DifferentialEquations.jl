
<a id='Extra-Functions-1'></a>

# Extra Functions

<a id='DifferentialEquations.modulechildren' href='#DifferentialEquations.modulechildren'>#</a>
**`DifferentialEquations.modulechildren`** &mdash; *Function*.



`modulechildren(m::Module)`

Returns the modules in m

<a id='DifferentialEquations.checkIfLoaded' href='#DifferentialEquations.checkIfLoaded'>#</a>
**`DifferentialEquations.checkIfLoaded`** &mdash; *Function*.



`checkIfLoaded(pkg::AbstractString)`

Returns true if module "pkg" is defined in Main, otherwise false.

<a id='DifferentialEquations.getNoise' href='#DifferentialEquations.getNoise'>#</a>
**`DifferentialEquations.getNoise`** &mdash; *Function*.



getNoise(N,node,elem;noiseType="White")

Returns a random vector corresponding to the noise type which was chosen.

<a id='DifferentialEquations.numparameters' href='#DifferentialEquations.numparameters'>#</a>
**`DifferentialEquations.numparameters`** &mdash; *Function*.



numparameters(f)

Returns the number of parameters of `f` for the method which has the most parameters.

