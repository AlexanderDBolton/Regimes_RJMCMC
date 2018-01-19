### The PELT Algorithm

The Pruned Exact Linear Time (PELT) algorithm was created by R. Killick, P. Fearnhead and I.A. Eckley. Their paper on the algorithm can be found at https://arxiv.org/pdf/1101.1438.pdf. Given data and a model to fit to the data (e.g. Bernoulli process, Markov chain, Poisson process, Gaussian) PELT finds the set of change points that minimises a given cost function. In this implementation the available cost functions are AIC and BIC. In the worst case the algorithm takes O(n^2) time (where n is the number of data points) but the algorithm will take O(n) time if the change points occur linearly with n.


### Input parameters

The input parameters must be stored in a file containing each of the following (see the Example directory for an example).

````
data_file
criterion
end_time
diff
````

- A data file: Should have a name of the form "yourPrefix_everything.txt". As for the RJMCMC sampler, the last 15 characters from this name will be removed and the remaining string will be used as a prefix to output files. The data file should give each of the simultaneous random processes listed on a separate line, e.g. sample_everything.txt could contain the following 4 lines:

````
4     markov_chain      sample_data_process_0.txt     uniform_arrivals
2     independent_chain sample_data_process_1.txt     sample_data_times_process_1.txt
3     gaussian    sample_data_process_2.txt     uniform_arrivals
3     poisson     sample_data_process_3.txt     sample_data_times_process_3.txt
````

The first number is the number of states: for a Markov chain this is the number of states (assuming that the states are labelled 0, 1, ...) in the input data. For a multinomial process (labelled independent_chain here) it is also the number of states (again assuming that the states are labelled from 0). For Gaussian observations 3 is used because the number of observations, sum of observations and sum of observations squared are the 3 sufficient statistics extracted from the data in any likelihood calculation. For Poisson observations 3 is used because, the number of observations, the sum of observations and the factorials of the observations are stored in order to calculate the log likelihood.

The second column gives the name of the process. "markov_chain", "independent_chain", "gaussian" or "poisson" can be used.

The third column gives the observations. They should be stored in separate files. The data values can be separated by tabs, spaces or put on different lines.

The third column gives a data file containing the times at which the data are observed. The observations must arrive in discrete time. Multiple observations can occur at the same time. If you supply your own file with observation times then they must be sorted in ascending order. The times at which observations can arrive are 0, 1, 2, ..., end_time (end_time is one of the parameters that you will specify). If your observations arrive regularly at times 0, 1, 2, ... then you can give "uniform_arrivals" as the input here and you don't need to create a file with observation times.
- criterion: This can be "AIC" or "BIC" and gives the information criterion that will be minimised. The set of change points that minimises the chosen information criterion will be returned.
- end_time: The time of the last arrival. Remember that the first observation is at time 0, so if you have n observations at times 0, 1, 2, ... then end_time should be n - 1.
- diff: The possible change point positions will change from 1, 2, ..., end_time to diff, 2*diff, .... Using a diff > 1 will reduce the amount of memory required to store the data and may improve the performance of the sampler if the random process is particularly long. The default value is 1.


### Output

The output will be a set of change points in a file called "yourPrefix_cp_n.txt".
