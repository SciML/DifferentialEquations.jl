# Notes on Algorithms

This page is a supplemental page which details some facts about the chosen algorithms,
why some I took the time to make optimized versions for, and for others why they
were ignored.

## Explicit Runge-Kutta ODE Algorithms

From what I can tell, this is by far the most comprehensive comparison of
Explicit Runge-Kutta ODE algorithms that you'll find.

## Order 4-

At this stage, coefficient of the truncation error seems to win out, or you are
willing to live with low tolerance anyways. Thus Bogacki-Shampine is the clear
winner in this category because at order 2/3 with FASL it has minimal numbers
of function evaluations but also is stable enough to step as needed. All other
methods don't compare because of the FASL property boosting the order and thus
the stability (for low orders, it pretty much holds that higher order = higher
stability (for optimal number of steps), which is not true as we go higher), making
it more stable and have less error for lower numbers of function evaluations than
the others in this category.

### Order 5

[Note that for all of these Peter Stone's modifications do not seem to be helpful since,
although they lower the truncation error, they also modify the stability region
in ways that can be worrisome (mostly they shrink the stability in the complex axis
near the origin, making the problems not as suitable for a "general purpose default"
like one would hope with a 4/5 solver)]

The "clear choice" is the Dormand-Prince 4/5 pair. This is the pair which is
used by default as ode45 in MATLAB, and serves similar functions in scipy,
ODE.jl, etc. The golden standard implementation is Hairer's DOPRI5 (offered
by ODEInterface.jl). After optimizations, DifferentialEquations.jl's native
DP5 solver is much more efficient (between 4x-400x) than DOPRI5's, with various design choices
factoring into this (which are documented in the benchmarks). This is pre-threading,
and within method threading will likely be at least doubled or tripled when
threading is enabled. Thus it's clear that the reference implementation to
try other methods against is the DifferentialEquations.jl DP5 method.

It's obvious that anything before Dormand-Prince 4/5's pair is simply not as
good because of the optimizations on the local truncation error coefficient
and the fact that FASL schemes essentially have one less function evaluation.
So the previous algorithms were implemented as tableaus for the historical
reasons but dealt with no further. These methods include the Runge, Cassity,
Butcher, Fehlburg, Lawson, Luther and Konen, and Kutta schemes.

The next set of schemes are the Papakostas-Papageorgiou schemes. The problem is
that they don't really get the much lower on the error than DP5, but also have
wacky stability near the origin.

Tsitouras's looks to be a good match against DP5 as a 6-stage scheme to take on
DP5. Its stability is similar to DP5 but its first error term is an order of magnitude
smaller. More tests will likely determine that this is much better than DP5 in
accordance with his paper.

Lastly, there are the 7-stage schemes. The more recent one is due to Sharp and Smart,
but I am ignoring this because its error term is almost an order of magnitude
larger than the BS pair, and its stability reagion is wonky near the origin.
Its only plus over the BS pair is that it has a slightly larger stability in the
real axis, which is not important when paired with adaptive stepping and for
use on non-stiff problems.

That leaves us with the Bogacki-Shampine pair. This pair gets more than an order
of magnitude lower truncation error, enhanced complex stability, and two error estimators
to make it more robust. In fact, this is the default which is chosen in
Mathematica. Its downside is that since it is an 8-stage scheme, it requires
an additional function evaluation.

Further tests will likely narrow this down to Tsitouras vs Bogacki-Shampine.
Who will come out on top? Who knows.

### Order 6

Sharp-Verner has bad complex stability near the origin. I don't like any of the
Peter Stone modifications here. Butcher and Chummund methods have stability issues near the
origin as well. Huta's method has too high of an error coefficient. Verner's 1991
has bad complex stability. Same as the most robust. The Verner "most efficient"
has really good stability and error coefficient. In fact, nothing is even
close except for Tsitouras' method. The DP method is two orders of magnitude
higher in error coefficient than Verner. The Luther methods have too much error.
Same as Tsitouras-Papakostas and  M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki.

Without a doubt the winner is the Verner "most efficient".

### Order 7

The Enright-Verner and other Verner methods all have stability issues near the
origin in the complex plane and higher error coefficients. Sharp and Smart
have higher error coefficients. Peter Stone's methods all have higher error.
It's very clear that the best here is the Tanaka-Yamashita (efficient, not the stable)
method by far.

### Order 8

The Cooper-Verner methods do not have an error estimate and so no adaptive
timestepping can be done. This is a deal-breaker. Going into this one would
think that the clear winner would be Dormand-Prince 8. But that's not the case.
However, that's comparing the classical 1981 DP87. Notice that the code for
Dop853 is based off of the 1989 paper which has different coefficients (and currently
I have no analysis for this).

The other methods include Verner's Maple dverk78 which is bested in both stability
and error coefficient by Enright-Verner which is bested by Tsitouras-Papakostas.

Thus the final showdown is between DP853 vs the Tsitouras-Papakostas pair.

### Order 9

The Tsitouras scheme and the Sharp scheme have funky stability near the origin.
Verner's schemes are much safer, and with similar error. They clearly
dominate this category.

### Order 10

Curtis' scheme has more function evalulations than needed, and Peter Stone's
modification reduces the truncation error by a lot but adds three more function
evaluations. Thus Hairer's 17 stage scheme (whose error and stability is similar
to Curtis') is clearly better. Once again Peter Stone's modification adds three
steps but does not reduce the truncation error here, so the unmodified version
does better.

Tom Baker's method increases the stability region to something which is more than
necessary but adds 4 function evaluations to do so (without lowering the error
very much). Ono's scheme minimizes the error more than Hairer's here, with
all else being basically the same. The Peter Stone methods add a lot of function
evalulations (5+) and so they would only be useful in the case where the function
evaluations are quick yet you still want really small error. Even then I'm not
convinced they are better than the other methods, or better than the higher order
methods which use less steps. The stability is only okay.

The Feagin scheme is fine, but with more error and less stability than the Hairer
scheme. Thus it seems clear that Hairer's method dominates this category.

### Order 11

The order 11 schemes are due to Tom Baker at the University of Teeside. They
have a nice sparsity pattern and receive slightly lower truncation error coefficents
than the Feagin, but Feagin's dominates by being "almost order 13" anyways
so while a nice try the order 11 scheme is likely overwhelmed in any case where
it would be useful.

### Order 12

Here there are the Feagin schemes and Ono's scheme. Ono's scheme gets horrible
stability with more error and so it's not in the running. Peter Stone's modifications
do not make a substantive change, and where they do they get rid of the nice
property that the Feagin 12 method satisfies many of the higher order conditions
as well, making it look even higher order on some problems. Thus the standard
Feagin 12 seems to win out in this category.

### Order 14

In this category there is just the Feagin. Peter Stone's modification barely
changes anything in the analysis so I did not even attempt it.
