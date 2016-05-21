using DifferentialEquations

SRIW1 = constructSRIW1()
bool1 = minimum(checkSRIOrder(SRIW1))

SRA1 = constructSRA1()
bool2 = minimum(checkSRAOrder(SRA1))

bool1 && bool2
