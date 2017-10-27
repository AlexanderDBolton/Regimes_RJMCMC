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

- A data file: 



notes gtools R package
      move .py file, PELT, RJMCMC, etc from hustler to git, create a new R file for analysing the results and get python to call it.
