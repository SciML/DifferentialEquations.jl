using DifferentialEquations

function test_num_parameters()
  f(x) = 2x
  f(x,y) = 2xy

  numpar = numparameters(f) # Should be 2

  numpar == 2
end
test_num_parameters()
