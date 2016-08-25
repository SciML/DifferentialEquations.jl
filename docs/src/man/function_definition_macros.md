# Function Definition Macros

DifferentialEquations.jl provides a set of macros for more easily and legibly
defining your differential equations. It exploits the standard notation for
mathematically writing differential equations and the notation for "punching
differential equations into the computer"; effectively doing the translation
step for you. This is best shown by an example. Say we want to solve the
[ROBER model](http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Single/DemoRobertson/demorobertson.pdf).
Using the `@ode_define` macro, we can do this by writing:

```julia
f = @ode_define begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁=>0.04 k₂=>3e7 k₃=>1e4
```

This looks just like psudocode! The macro will expand this to the "standard form",
i.e. the ugly computer form:

```julia
f = (t,u,du) -> begin
  du[1] = -0.04*u[1] + 1e4*u[2]*u[3]
  du[2] = 0.04*u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3]
  du[3] = 3e7*u[2]^2
end
```

The other macro which is currently provided is the `@fem_define` macro. This macro
is for parsing and writing FEM functions. For example, in the FEM methods you have
to use `x[:,1]` instead of `x` and `x[:,2]` instead of `y`. The macro will automatically
do this replacement, along with adding in parameters. Since FEM functions are more
general, we also have to give it the function signature. Using the macro looks like this:

```julia
f  = @fem_define((x),(),begin
  sin(α.*x).*cos(α.*y)
end,α=>2π)
gD = @fem_define((x),(),begin
  sin(α.*x).*cos(α.*y)/β
end,α=>2π,β=>8π*π)
```

This is equivalent to the definition

```julia
f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
gD(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
```

The true power comes in when dealing with nonlinear equations. The second argument,
which we skipped over as `()`, is for listing the variables you wish to define the
equation by. Mathematically you may be using `u`,`v`,`w`, etc., but for array-based
algorithms you need to use `u[:,1]`,`u[:,2]`,etc. To avoid obfuscated code, the
`@fem_define` macro does this conversion. For example:

```julia
l = @fem_define((t,x,u),(u,v),begin
  [ones(length(u))-α*u ones(length(v))-v]
end,α=>0.5)
```
says there are two equations, one for `u` (`ones(length(u))-α*u`) and one for `v`
`(ones(length(v))-v)`. This expands to the equation

```julia
h = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]   ones(size(x,1))-u[:,2]]
```

When you have 10+ variables, using `@fem_define` leads to code which is much
easier to read!
