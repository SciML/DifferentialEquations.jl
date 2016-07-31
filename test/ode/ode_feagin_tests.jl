using DifferentialEquations
srand(100)
u = [0.5,0.2]

a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(eltype(u))

typeof(a0100) == Float64
eltype(b) == Float64
eltype(c) == Float64
u2 = [BigFloat(0.5),BigFloat(0.2)]

a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(eltype(u2))
typeof(a0100) == BigFloat
eltype(b) == BigFloat
eltype(c) == BigFloat

prob = twoDimlinearODEExample()

## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(12:-1:4)
testTol = 0.2

println("Feagin RKs")
sol =solve(prob::ODEProblem,Δt=Δts[1],alg=:Feagin10)

prob = twoDimlinearODEExample(α=ones(BigFloat,4,2),u₀=map(BigFloat,rand(4,2)).*ones(4,2)/2)

sim = test_convergence(Δts,prob,alg=:Feagin10)
plot(sim); Plots.gui()
#sim = test_convergence(Δts,prob,alg=:RK4)
bool1 = abs(sim.𝒪est[:final]-10) < testTol
