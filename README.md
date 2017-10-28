# Regimes_RJMCMC
C++ Code for performing unmarked, binary and regime change points RJMCMC

Language Requirements

- C++ (with the Gnu Scientific Library (GSL) installed),
- Python (to recreate Figure 8 in Malware Family Discovery Using Reversible Jump MCMC Sampling of Regimes), and
- R (with the gtools package installed in order to recreate Figure 8).


Code

The code for performing RJMCMC sampling is written in C++ and is contained in the 'Code' directory. When run in terminal on a computer running Linux, the compile_regime executable will compile the code into an executable called RJMCMC. See below (the sections titled 'Input Parameters' and 'How To Run The RJMCMC Sampler') for information on how to run a job.


Example

The example gives two simultaneous Bernoulli processes and a full description of the change points for these processes. Input parameters are included as well.


Simulation Study

The Simulation_Study folder contains two executables: RJMCMC and PELT, which run the algorithms after which they are named. Running the .py file will generate the 1000 samples described in Section 6.2 of Malware Family Discovery Using Reversible Jump MCMC Sampling of Regimes. By default the generation of samples and the running of the algorithm will run on 4 cores, but this can be changed in the call to "Pool(processes = 4)".


Input Parameters

The input parameters must be stored in a file containing each of the following (see the Example directory for an example).
- A data file: Should have a name of the form "yourPrefix_everything.txt" as the last 15 characters from this name will be removed and the remaining string will be used as a prefix to output files. The data file should detail each of the simultaneous random processes listed on a separate line, e.g. sample_everything.txt could contain the following 4 lines:
4     markov_chain      sample_data_process_0.txt     uniform_arrivals
2     independent_chain sample_data_process_1.txt     sample_data_times_process_1.txt
3     gaussian    sample_data_process_2.txt     uniform_arrivals
2     poisson     sample_data_process_3.txt     sample_data_times_process_3.txt
The first number is the number of states: for a Markov chain this is the number of states (assuming that the states are labelled 0, 1, ...) in the input data. For a multinomial process (labelled independent_chain here) it is also the number of states (again assuming that the states are labelled from 0). For Gaussian observations 3 is used because the number of observations, sum of observations and sum of observations squared are the 3 sufficient statistics extracted from the data in any likelihood calculation. For Poisson observations 2 is used because only the number of observations and the sum of observations are used in any likelihood calculation.
The second column gives the name of the process. markov_chain, independent_chain, gaussian or poisson can be used.
The third column gives the observations. They should be stored in separate files. The data values can be separated by tabs, spaces or put on different lines.
The third column gives a data file containing the times at which the data are observed. The observations must arrive in discrete time. Multiple observations can occur at the same time. If you supply your own file with observation times then they must be sorted in ascending order. The times at which observations can arrive are 0, 1, 2, ..., end_time (end_time is one of the parameters that you will specify). If your observations arrive regularly at times 0, 1, 2, ... then you can give "uniform_arrivals" as the input here and you don't need ot create a file with observation times.
- end_time: The time of the last arrival. Remember that the first observation is at time 0, so if you have n observations at times 0, 1, 2, ... then end_time should be n - 1.
- diff: The possible change point positions will change from 1, 2, ..., end_time to diff, 2*diff, .... Using a diff > 1 will reduce the amount of memory required to store the data and may improve the performance of the sampler if the random process is particularly long. The default value is 1.
- gaussian_mu_0: The notation for the Gaussian parameters is taken from 'Conjugate Bayesian analysis of the Gaussian distribution' by Kevin P. Murphy http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf The prior is Normal-Gamma, so the joint prior for the mean mu and precision lambda (the precision is the inverse of the variance) is mu ~ N(mu_0, 1/(kappa_0 * lambda)), lambda ~ Gamma(alpha_0, rate = beta_0). An uninformative value for gaussian_mu_0 is 0.
- gaussian_kappa_0: Low values of kappa_0 are uninformative.
- gaussian_alpha_0: Low values of alpha_0 are uninformative.
- gaussian_beta_0: Low values of beta_0 are uninformative.


notes gtools R package
      move .py file, PELT, RJMCMC, etc from hustler to git, create a new R file for analysing the results and get python to call it.
