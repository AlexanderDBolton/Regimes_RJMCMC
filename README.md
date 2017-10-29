# Regimes_RJMCMC
## C++ Code for performing unmarked, binary and regime change points RJMCMC

### Language Requirements

- C++ with the Gnu Scientific Library (GSL) installed,
- Python (if you want to recreate Figure 8 in Malware Family Discovery Using Reversible Jump MCMC Sampling of Regimes), and
- R with the gtools package installed (if you want to recreate Figure 8).


### Code

The code for performing RJMCMC sampling is written in C++ and is contained in the 'Code' directory. When run in terminal on a computer running Linux, the compile_regime executable will compile the code into an executable called RJMCMC. See below (the sections titled 'Input Parameters' and 'How To Run The RJMCMC Sampler') for information on how to run a job.


### Example

The example gives two simultaneous Bernoulli processes and a full description of the change points for these processes. Input parameters are included as well.


### Simulation Study

The Simulation_Study folder contains two executables: RJMCMC and PELT, which run the algorithms after which they are named. Running the .py file will generate the 1000 samples described in Section 6.2 of Malware Family Discovery Using Reversible Jump MCMC Sampling of Regimes. By default the generation of samples and the running of the algorithm will run on 4 cores, but this can be changed in the call to "Pool(processes = 4)".


### Input Parameters

The input parameters must be stored in a file containing each of the following (see the Example directory for an example).
- A data file: Should have a name of the form "yourPrefix_everything.txt" as the last 15 characters from this name will be removed and the remaining string will be used as a prefix to output files. The data file should detail each of the simultaneous random processes listed on a separate line, e.g. sample_everything.txt could contain the following 4 lines:

````
4     markov_chain      sample_data_process_0.txt     uniform_arrivals
2     independent_chain sample_data_process_1.txt     sample_data_times_process_1.txt
3     gaussian    sample_data_process_2.txt     uniform_arrivals
2     poisson     sample_data_process_3.txt     sample_data_times_process_3.txt
````

The first number is the number of states: for a Markov chain this is the number of states (assuming that the states are labelled 0, 1, ...) in the input data. For a multinomial process (labelled independent_chain here) it is also the number of states (again assuming that the states are labelled from 0). For Gaussian observations 3 is used because the number of observations, sum of observations and sum of observations squared are the 3 sufficient statistics extracted from the data in any likelihood calculation. For Poisson observations 2 is used because only the number of observations and the sum of observations are used in any likelihood calculation.

The second column gives the name of the process. "markov_chain", "independent_chain", "gaussian" or "poisson" can be used.

The third column gives the observations. They should be stored in separate files. The data values can be separated by tabs, spaces or put on different lines.

The third column gives a data file containing the times at which the data are observed. The observations must arrive in discrete time. Multiple observations can occur at the same time. If you supply your own file with observation times then they must be sorted in ascending order. The times at which observations can arrive are 0, 1, 2, ..., end_time (end_time is one of the parameters that you will specify). If your observations arrive regularly at times 0, 1, 2, ... then you can give "uniform_arrivals" as the input here and you don't need ot create a file with observation times.
- end_time: The time of the last arrival. Remember that the first observation is at time 0, so if you have n observations at times 0, 1, 2, ... then end_time should be n - 1.
- diff: The possible change point positions will change from 1, 2, ..., end_time to diff, 2*diff, .... Using a diff > 1 will reduce the amount of memory required to store the data and may improve the performance of the sampler if the random process is particularly long. The default value is 1.
- gaussian_mu_0: The notation for the Gaussian parameters is taken from 'Conjugate Bayesian analysis of the Gaussian distribution' by Kevin P. Murphy http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf The prior is Normal-Gamma, so the joint prior for the mean mu and precision lambda (the precision is the inverse of the variance) is mu ~ N(mu_0, 1/(kappa_0 * lambda)), lambda ~ Gamma(alpha_0, rate = beta_0). An uninformative value for gaussian_mu_0 is 0. If no processes with Gaussian data are run, then this can be set to be any value.
- gaussian_kappa_0: Low values of kappa_0 are uninformative. If no processes with Gaussian data are run, then this can be set to be any value.
- gaussian_alpha_0: Low values of alpha_0 are uninformative. If no processes with Gaussian data are run, then this can be set to be any value.
- gaussian_beta_0: Low values of beta_0 are uninformative. If no processes with Gaussian data are run, then this can be set to be any value.
- poisson_alpha: The prior for the Poisson parameter lambda is Gamma(poisson_alpha, poisson_beta). If no processes with Poisson data are run, then this can be set to be any value.
- poisson_beta: If no processes with Gaussian data are run, then this can be set to be any value.
- markov_alpha: For a multinomial process or a Markov chain, the prior distribution is Dirichlet(markov_alpha, ..., markov_alpha). If no multinomial processes or Markov chains are run, then this can be set to be any value.
- p: The expected number of change points. This will be converted to a prior of 'p = expected number of change points / possible number of change points' and each possible change point position has prior probability p of being a change point.
- var_p: If var_p is set to be 0 then the prior distribution of the change points is a Bernoulli process, with each point having fixed probability p of being a change point. If 0 < var_p < p*(1 - p) then the prior is a beta distribution with mean p and variance var_p (the beta distribution parameters are estimated using the method of moments.
- beta_alpha: The prior for the probability that a change point affects a process is Beta(beta_alpha, beta_alpha).
- dirichlet_alpha: The prior for which regime begins at a change point is Dirichlet(dirichlet_alpha, ..., dirichlet_alpha).
- rho: The prior distribution of the number of regimes is Geometric(rho).
- burnin: The number of burn-in iterations to run before recording simulations at each of the three stages.
- iterations: The number of RJMCMC iterations to be run in each of the three stages.
- thinning: The number of iterations to skip between iterations. The total number of iterations recorded is iterations / thinning.
- number_of_association_matrix_bins: To build an association matrix (a matrix that states the posterior probability that a pair of time points in a process are generated from the same regime), the number of bins that the points should be divided into.
- record_basic_rj: Input 0 or 1. If 0, the output of the unmarked change point RJMCMC sampling will not be recorded. If 1, it will be recorded.
- record_binary_rj: Input 0 or 1. If 0, the output of the binary marked vector and change point RJMCMC sampling will not be recorded. If 1, it will be recorded.
- record_full_rj: Input 0 or 1. If 0, the output of the regime marked vector and change point RJMCMC sampling will not be recorded. If 1, it will be recorded.
- record_similarity_matrix: Input 0 or 1. If 0, the matrix giving the similarities of the different traces will not be included. If 1, it will be included.
- starting_changepoints_file: Gives the name of the file containing a starting set of change points. If no set of starting change points is provided, write "no_starting_changepoints".
- separator_file: The name of a file containing the times at which the different traces begin (0 should not be included). If no separators are included then write "no_separators".
- trace_lengths_file: The name of a file containing the number of observations in each trace. Use if diff > 1 as this will make the trace length estimates incorrect. Write "no_trace_lengths" if diff == 1.
- seed: The random seed for the RJMCMC sampler.


### How to Run The RJMCMC Sampler
- Create a yourPrefix_input_parameters.txt file using the parameters above along with yourPrefix_everything.txt and data files.
- Assuming that the executable RJMCMC is in the same directory as your data and input parameters, run
````
./RJMCMC yourPrefix_input_parameters.txt
````

notes gtools R package
      remove "Assuming... " line from basic_RJMCMC.h and then recompile, move executable back into GitHub. move .py file, PELT, RJMCMC, etc from hustler to git, create a new R file for analysing the results and get python to call it.
