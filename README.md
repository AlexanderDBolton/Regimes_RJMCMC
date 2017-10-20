# Regimes_RJMCMC
C++ Code for performing unmarked, binary and regime change points RJMCMC

Language Requirements

- C++ (with the Gnu Scientific Library (GSL) installed),
- Python (to recreate Figure 8 in Malware Family Discovery Using Reversible Jump MCMC Sampling of Regimes), and
- R (with the gtools package installed in order to recreate Figure 8).


Code

The code for performing RJMCMC sampling is written in C++ and is contained in the 'Code' directory. The compile_regime executable will compile the code into an executable called RJMCMC.




Simulation Study

The Simulation_Study folder contains two executables: RJMCMC and PELT. These run the algorithms after which they are named. 

notes gtools R package
      move .py file, PELT, RJMCMC, etc from hustler to git, create a new R file for analysing the results and get python to call it.
