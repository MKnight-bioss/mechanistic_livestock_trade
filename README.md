# mechanistic_livestock_trade
C++ code of model described in "A mechanistic model captures livestock trading, disease dynamics, and compensatory behaviour in response to control measures"

Code requires installation of GSL random number generator. This can be obtained from https://www.gnu.org/software/gsl/ or code can be edited to use user desired random number generator. Users can replace any calls to the gsl rng with their desired rng by removing code where necessary. When compiling, add `-lgsl -lgslcblas`

Code requires a single input file, `trade_quantities_final.csv`. This file contains stock generation rates, average trade quantities, and values for the a, b, and partnership cessation rate required to simulate the system. All code is presented within the `main` function. Users can freely edit code as desired.

The C++ file `trade_model_demand_shocks.cpp` performs the simulations analysing the system response to shocks in farm-level demand. This will produce two output files, one outputting time series data of trade quantities for each shock, the other outputting system-average values of trade quantities for each shock.

The C++ file `trade_model_changes_to_epsilon_b.cpp` performs the simulations analysing the system with changes to trading propensities via different values of epsilon_b as described in the main text. Running this model will produce two output files, one outputting time series data of trade quantities and disease prevalence for each value of epsilon_b, the other outputting equilibrium system-average values of trade quantities and disease prevalence for each value of epsilon_b.

The C++ file `trade_model_whole_batch_testing.cpp` performs the simulations analysing the system in response to the introduction of whole batch testing and rejection if an animal in the traded batch tests positive for infection. Two output files are produced. One outputting time series data of the trading system and disease prevalence for various values of the test sensitivity, the other outputting the long-run equilibrium of the trading system and disease prevalence for the same values of the test sensitivity.
