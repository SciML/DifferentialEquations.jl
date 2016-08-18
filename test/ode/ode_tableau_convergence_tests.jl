# This definitely needs cleaning
using DifferentialEquations
probArr = Vector{DEProblem}(2)
bigprobArr = Vector{DEProblem}(2)
probArr[1] = linearODEExample()
probArr[2] = twoDimlinearODEExample()
bigprobArr[1] = linearODEExample(u₀=BigFloat(1),α=BigFloat(1))
bigprobArr[2] = twoDimlinearODEExample(α=ones(BigFloat,4,2),u₀=map(BigFloat,rand(4,2)).*ones(4,2)/2)
srand(100)
## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(8:-1:4)
testTol = 0.2
superduperbool = Vector{Bool}(2)
bools = Vector{Bool}(0)
i=1
#for i = 1:2
prob = probArr[i]
bigprob = bigprobArr[i]


# Order 5
Δts = 1.//2.^(7:-1:4)
tab = constructRungeFirst5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructCassity5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructLawson5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructLutherKonen5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructLutherKonen52()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructLutherKonen53()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructPapakostasPapaGeorgiou5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructPapakostasPapaGeorgiou52()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructTsitouras5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

Δts = 1.//2.^(6:-1:4)
tab = constructBogakiShampine5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

tab = constructSharpSmart5()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

# Order 6

tab = constructButcher6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructButcher62()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructButcher63()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

Δts = 1.//2.^(5:-1:2)
tab = constructDormandPrince6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructSharpVerner6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructVerner916()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructVerner9162()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructVernerRobust6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructVernerEfficient6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructPapakostas6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructLawson6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructTsitourasPapakostas6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructDormandLockyerMcCorriganPrince6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructTanakaKasugaYamashitaYazaki6D()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructTanakaKasugaYamashitaYazaki6C()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructTanakaKasugaYamashitaYazaki6B()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructTanakaKasugaYamashitaYazaki6A()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructMikkawyEisa()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructChummund6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructChummund62()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructHuta6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructHuta62()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructVerner6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructDverk()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

tab = constructClassicVerner6()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

# Order 7

tab = constructButcher7()
sim = test_convergence(Δts,bigprob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)

Δts = 1.//2.^(5:-1:2)
tab = constructClassicVerner7()
sim = test_convergence(Δts,bigprob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)

# Order 8
Δts = 1.//2.^(4:-1:1)
tab = constructClassicVerner8()
sim = test_convergence(Δts,bigprob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)

Δts = BigFloat(1.)./BigFloat(2.0).^(10:-1:7)
# Order 9

tab = constructClassicVerner9()
sim = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)

#end
