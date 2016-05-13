# Overview of DifferentialEquations.jl Usage

The general workflow for using the package is as follows:

* Define a problem
* Generate a mesh
* Use a solver on the problem and mesh
* Analyze the output

Problems are specified via a type interface. For example, for the Poisson equation
``Δu = f``, one defines a type which holds `f` and the boundary condition functions.

Next, one generates a mesh. For example, if one wants to solve the Heat equation
in the parabolic cylinder of the unit square, i.e. [0,1]^2 x [0,T], then one
has to discretize this. Tools within the package will generate meshes from
general characteristics. For example, most tools require only specifying the
general shape, Δx, Δt, and T and will generate the mesh.

One then passes the mesh and the problem to the solver interface. The solver then
solves the differential equation using the some numerical methods (which can be
specified via keyword arguments). The solver returns a solution object which
hold all of the details for the solution.

With the solution object, you do the analysis as you please! For some result `res`,
the field `res.u` returns the final solution, and if you give a true solution,
`res.uTrue` is the true solution at the final time. If you specified to the solver
`fullSave=true`, then `res.uFull` and `res.tFull` will be outputted which hold the
solution/time at every `saveSteps` (default set to 100, meaning it saves an output
every 100 steps).

However, DifferentialEquations.jl also provides some helper functionality to assist
with general forms of analysis. Various plotter functions take in result types
and use these to make standard plots. If `fullSave` was used, the plotters can
generate animations of the solutions to evolution equations.  An array of solutions
can be made into a `ConvergenceSimulation` which then generates all of the
convergence test results and allows for plotting (great for developing new methods!).
