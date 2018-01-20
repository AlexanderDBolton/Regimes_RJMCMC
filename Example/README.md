The data in this directory are two Bernoulli processes subject to changes in regimes at times 5,000, 10,000 and 15,000. One process has some censored values so not all observations are available (and the censoring occurs completely at random).

To run the RJMCMC sampler on the example in this directory in Linux, run
````
./RJMCMC example_input_parameters.txt
````
and to run the PELT algorithm to detect the set of change points that minimises the BIC, run
````
./PELT example_PELT_input_parameters.txt
````
