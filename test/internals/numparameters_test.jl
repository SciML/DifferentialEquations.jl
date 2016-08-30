using DifferentialEquations

function test_num_parameters()
  f(x) = 2x
  f(x,y) = 2xy

  numpar = numparameters(f) # Should be 2
  g = (x,y) -> x^2
  numpar2 = numparameters(g)
  numpar == 2 && numpar2 == 2
end
test_num_parameters()
