
<a id='The-Solution-Type-1'></a>

# The Solution Type


appxTrue!


<a id='Related-Functions-1'></a>

## Related Functions

<a id='DifferentialEquations.FEMSolution' href='#DifferentialEquations.FEMSolution'>#</a>
**`DifferentialEquations.FEMSolution`** &mdash; *Type*.



FEMSolution

A type which holds the data for the solution to a finite element problem.

<a id='DifferentialEquations.appxTrue!' href='#DifferentialEquations.appxTrue!'>#</a>
**`DifferentialEquations.appxTrue!`** &mdash; *Function*.



appxTrue!(res,res2)

Adds the solution from res2 to the FEMSolution object res. Useful to add a quasi-true solution when none is known by computing once at a very small time/space step and taking that solution as the "true" solution

<a id='DifferentialEquations.PdeSolution' href='#DifferentialEquations.PdeSolution'>#</a>
**`DifferentialEquations.PdeSolution`** &mdash; *Type*.



PdeSolution: Wrapper for the objects obtained from a PdeSolver

