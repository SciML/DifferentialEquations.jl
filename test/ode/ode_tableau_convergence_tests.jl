# This definitely needs cleaning
using DifferentialEquations
probArr = Vector{DEProblem}(2)
bigprobArr = Vector{DEProblem}(2)
probArr[1] = linearODEExample()
probArr[2] = twoDimlinearODEExample()
bigprobArr[1] = linearODEExample(u₀=BigFloat(1),α=BigFloat(1))
bigprobArr[2] = twoDimlinearODEExample(α=ones(BigFloat,4,2),u₀=map(BigFloat,rand(4,2)).*ones(4,2)/2)
setprecision(400)
srand(100)
## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(8:-1:4)
testTol = 0.3
superduperbool = Vector{Bool}(3)
bools = Vector{Bool}(0)

for i = 1:3
  if i>1
    prob = probArr[2]
    bigprob = bigprobArr[2]
  else
    prob = probArr[1]
    bigprob = bigprobArr[1]
  end
  if i>2
    alg = :ExplicitRKVectorized
  else
    alg = :ExplicitRK
  end

  # Order 2

  tab = constructHuen()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-2) < testTol)

  tab = constructRalston()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-2) < testTol)

  # Order 3

  tab = constructBogakiShampine3()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-3) < testTol)

  # Order 4



  # Order 5

  Δts = 1.//2.^(7:-1:4)
  tab = constructRKF5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  Δts = 1.//2.^(7:-1:4)
  tab = constructDormandPrince()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  sim2 = test_convergence(Δts,prob,alg=:DP5)
  push!(bools,(abs(sim.𝒪est[:l∞]-5) < testTol) && (maximum(sim[end][end]-sim2[end][end]) < 1e-10))

  tab = constructCashKarp()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  Δts = 1.//2.^(7:-1:4)
  tab = constructRungeFirst5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructCassity5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructLawson5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructLutherKonen5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructLutherKonen52()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructLutherKonen53()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructPapakostasPapaGeorgiou5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructPapakostasPapaGeorgiou52()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructTsitouras5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  Δts = 1.//2.^(6:-1:4)
  tab = constructBogakiShampine5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)

  tab = constructSharpSmart5()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5) < testTol)


  # Order 6

  Δts = 1.//2.^(6:-1:4)
  tab = constructButcher6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructButcher62()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructButcher63()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  Δts = 1.//2.^(5:-1:1)
  tab = constructDormandPrince6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6.7) < testTol) # Better on linear

  tab = constructSharpVerner6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructVerner916()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructVerner9162()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructVernerRobust6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  #=
  tab = constructVernerEfficient6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)
  =#

  tab = constructPapakostas6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructLawson6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  Δts = 1.//2.^(3:-1:1)
  tab = constructTsitourasPapakostas6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6.7) < testTol) # Better on linear

  Δts = 1.//2.^(5:-1:1)
  tab = constructDormandLockyerMcCorriganPrince6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  #=
  tab = constructTanakaKasugaYamashitaYazaki6D()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)


  tab = constructTanakaKasugaYamashitaYazaki6C()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructTanakaKasugaYamashitaYazaki6B()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructTanakaKasugaYamashitaYazaki6A()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  =#

  tab = constructMikkawyEisa()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6.53) < testTol) # Odd behavior

  tab = constructChummund6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  #=
  tab = constructChummund62()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)
  =#

  Δts = 1.//2.^(4:-1:1)
  tab = constructHuta6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-5.5) < testTol) # Low convergence, error noted in Stone notes

  Δts = 1.//2.^(5:-1:1)
  tab = constructHuta62()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructVerner6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6.7) < testTol) # Better on linear

  Δts = 1.//2.^(4:-1:1)
  tab = constructDverk()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  tab = constructClassicVerner6()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6) < testTol)

  # Order 7

  Δts = 1.//2.^(5:-1:1)
  tab = constructButcher7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)

  Δts = 1.//2.^(5:-1:2)
  tab = constructClassicVerner7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)

  #=
  tab = constructVernerRobust7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  Δts = 1.//2.^(5:-1:1)
  tab = constructEnrightVerner7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7.15) < testTol) # Better on linear

  #=
  tab = constructTanakaYamashitaStable7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  #=
  tab = constructTanakaYamashitaEfficient7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  #=
  tab = constructSharpSmart7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  Δts = 1.//2.^(4:-1:1)
  tab = constructSharpVerner7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-6.5) < testTol) # Coefficients aren't accurate enough, drop off error

  #=
  tab = constructVernerEfficient7()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  # Order 8
  Δts = 1.//2.^(4:-1:1)
  tab = constructClassicVerner8()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)

  Δts = 1.//2.^(3:-1:1)
  tab = constructCooperVerner8()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7.75) < testTol) #Coefficients not accurate enough

  Δts = 1.//2.^(4:-1:1)
  tab = constructCooperVerner82()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)

  #=
  tab = constructTsitourasPapakostas8()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)
  =#

  #=
  Δts = 1.//2.^(4:-1:1)
  tab = constructdverk78()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)
  =#

  #=
  Δts = 1.//2.^(4:-1:1)
  tab = constructEnrightVerner8()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)
  =#

  Δts = 1.//2.^(4:-1:1)
  tab = constructCurtis8()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)

  Δts = 1.//2.^(4:-1:1)
  tab = constructRKF8(BigFloat)
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)

  #=
  tab = constructDormandPrince8(BigFloat)
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8) < testTol)
  =#

  Δts = 1.//2.^(3:-1:1)
  tab = constructDormandPrince8_64bit(BigFloat)
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-8.4) < testTol)

  # Order 9

  #=
  tab = constructClassicVerner9()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  #=
  tab = constructVernerRobust9()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  #=
  tab = constructVernerEfficient9()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  Δts = 1.//2.^(4:-1:1)
  #=
  tab = constructSharp9()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#
  #=
  tab = constructTsitouras9()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  #=
  tab = constructTsitouras92()
  sim = test_convergence(Δts,prob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-7) < testTol)
  =#

  ## Order 10

  Δts = 1.//2.^(5:-1:1)
  tab = constructCurtis10()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-10) < testTol) #Coefficients too crude

  tab = constructOno10()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-10) < testTol) #Coefficients too crude

  Δts = 1.//2.^(5:-1:1)
  tab = constructFeagin10Tableau()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-10) < testTol)

  tab = constructCurtis10()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-10) < testTol)

  Δts = 1.//2.^(6:-1:1)
  tab = constructBaker10()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-10.8) < testTol)


  tab = constructHairer10()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-10.8) < testTol)

  ## Order 12

  Δts = 1.//2.^(6:-1:1)
  tab = constructFeagin12Tableau()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-12.6) < testTol)

  tab = constructOno12()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-11.6) < testTol)

  ## Order 14

  Δts = 1.//2.^(6:-1:1)
  tab = constructFeagin14Tableau()
  sim = test_convergence(Δts,bigprob,alg=alg,tableau=tab)
  push!(bools,abs(sim.𝒪est[:l∞]-15.5) < testTol)

  superduperbool[i] = minimum(bools)
end

minimum(superduperbool)
