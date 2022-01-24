# mechanistic_livestock_trade
C++ code of model described in "A mechanistic model captures livestock trading, disease dynamics, and compensatory behaviour in response to control measures"

Code requires installation of GSL random number generator. This can be obtained from https://www.gnu.org/software/gsl/ or code can be edited to use user desired random number generator. Users can replace any calls to the gsl rng with their desired rng by removing code where necessary. When compiling, add `-lgsl -lgslcblas`

Code requires a single input file, "trade_quantities_final.csv". This file contains stock generation rates, average trade quantities, and values for the a, b, and partnership cessation rate required to simulate the system. 
