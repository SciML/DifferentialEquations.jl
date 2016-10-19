# Differential Algebraic Equation (DAE) Example

This tutorial will introduce you to the functionality for solving DAEs. Other
introductions can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/examples).

In this example we will solve the equation

```math
f(t,u,du) = 0
```

where `f` is the a variant of the Roberts equation. This equations is actually of
the form

```math
\begin{align}
du = f(t,u)  
 0 = g(t,u)
 \end{align}
```

or is also known as a constrained differential equation where `g` is the constraint
equation. The Roberts model can be written in the form:

```math
\begin{align}
dy_1 &= -0.04y₁ + 10^4 y_2 y_3
dy_2 &= 0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2
1 &=  y_{1}  y_{2} + y_{3}
\end{align}
```

with initial conditions ``y_1(0) = 1``, ``y_2(0) = 0``, ``y_3(0) = 0``,
``dy_1 = - 0.04``, ``dy_2 = 0.04``, and ``dy_3 = 0.0``.

The workflow for DAEs is the same as for the other types of equations, where all
you need to know is how to define the problem. A DAEProblem is specified by defining
an in-place update `f(t,u,du,out)` which uses the values to mutate `out` as the
output. To makes this into a DAE, we move all of the variables to one side.
Thus we can define the function:

```julia
f = function (t,u,du,out)
  out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
  out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
  out[3] = u[1] + u[2] + u[3] - 1.0
end
```

with initial conditons

```julia
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
```

and make the DAEProblem:

```julia
prob = DAEProblem(f,u₀,du₀)
```

As with the other DifferentialEquations problems, the commands are then to solve
and plot:

```julia
tspan = [0;100000]
sol = solve(prob,tspan)
plot(sol)
```

which, despite how interesting the model looks, produces a relatively simple
output:

![IntroDAEPlot](../assets/introDAEPlot.png)
