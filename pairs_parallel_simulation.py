import multiprocessing as mp
from multiprocessing import Pool
import os
import os.path
import math

def run_simulation(x):
    new_filename = str(x)
    # Run R code to generate the data files and the true cp bins, k, log-posterior
    os.system("Rscript pairs_simulation_study.R " + new_filename)
    
    # Create input files for running PELT in order to get the change points that minimise the BIC.
    with open(new_filename + "_PELT_input_parameters.txt", "w") as PELT_input_parameters:
        PELT_input_parameters.write(new_filename + "_everything.txt\n")
        PELT_input_parameters.write("BIC\n")
        PELT_input_parameters.write("19999\n")
        PELT_input_parameters.write("1")
    
    with open(new_filename + "_everything.txt", "w") as everything:
        everything.write("4\tindependent_chain" + "\t" + new_filename + "_data_state_" + str(0) + ".txt\t" + new_filename + "_times_state_" + str(0) + ".txt")
        for state in range(1, 4):
            everything.write("\n" + "4\tindependent_chain" + "\t" + new_filename + "_data_state_" + str(state) + ".txt" + "\t" + new_filename + "_times_state_" + str(state) + ".txt")

    with open(new_filename + "_separator_changepoint.txt", "w") as sep:
        sep.write("10000")
    
    # Run the PELT code to get the BIC change points
    os.system("./PELT " + new_filename + "_PELT_input_parameters.txt")

    # If 100,000 is in the PELT set of change points it needs to be removed because it is a separator
    with open(new_filename + "_cp_n.txt", "r") as input:
        with open(new_filename + "_PELT_cp_n.txt", "wb") as output:
            for line in input:
                if line != "10000\n":
                    output.write(line)

    # Create input files for the RJMCMC system
    with open(new_filename + "_input_parameters.txt", "w") as inputparameters:
        inputparameters.write(new_filename + "_everything.txt\n")
        inputparameters.write("19999\n")
        inputparameters.write("1\n")
        inputparameters.write("0\n1\n1\n1\n")
        inputparameters.write("1\n1\n")
        inputparameters.write("1\n")
        inputparameters.write("1\n0\n")
        inputparameters.write("1\n1\n")
        inputparameters.write("0.2\n")
        inputparameters.write("5000000\n5000000\n")
        inputparameters.write("1000\n")
        inputparameters.write("0\n0\n0\n1\n1\n")
        inputparameters.write(new_filename + "_PELT_cp_n.txt\n")
        inputparameters.write(new_filename + "_separator_changepoint.txt\n")
        inputparameters.write("no_trace_lengths\n")
        inputparameters.write("0")

    # Run the RJMCMC sampler on the concatenated traces
    os.system("./RJMCMC " + new_filename + "_input_parameters.txt")

    # Read in actual number of full and effective CPs and the min proportion similarity
    with open(new_filename + "_actual_number_of_CPs.txt", "r") as actual_k_file:
        actual_k = int(actual_k_file.readline()) # find the actual k value

    with open(new_filename + "_actual_effective_number_of_CPs.txt", "r") as actual_effective_k_file:
        actual_effective_k = int(actual_effective_k_file.readline())

    with open(new_filename + "_actual_min_proportion_similarity.txt", "r") as actual_min_proportion_similarity_file:
        actual_min_proportion_similarity = float(actual_min_proportion_similarity_file.readline())

    with open(new_filename + "_actual_similarity.txt", "r") as actual_similarity_file:
        actual_similarity = float(actual_similarity_file.readline())

    # Read in the distribution of the full and effective CPs and the min proportion similarity file
    full_k_dist = []
    with open(new_filename + "_full_dimension_distribution.txt", "r") as k_file:
        for line in k_file:
            prob = float(line.split("\t")[1])
            full_k_dist.append(prob)

    full_effective_k_dist = []
    with open(new_filename + "_full_effective_dimension_distribution.txt", "r") as k_file:
        for line in k_file:
            prob = float(line.split("\t")[1])
            full_effective_k_dist.append(prob)

    min_proportion_trace = []
    with open(new_filename + "_min_proportion_similarity_matrices.txt", "r") as min_prop_trace_file:
        next(min_prop_trace_file) # skip the first line, which is a header
        for line in min_prop_trace_file:
            sim = float(line)
            min_proportion_trace.append(sim)

    similarity_trace = []
    with open(new_filename + "_similarity_matrices.txt", "r") as similarity_trace_file:
        next(similarity_trace_file) # skip the first line, which is a header
        for line in similarity_trace_file:
            sim = float(line)
            similarity_trace.append(sim)

    # Compare log-post, k, binning
    with open("comparison_results_" + new_filename + ".txt", "w") as comparison_file:
        # Calculate the expected number of changepoints. Write the difference between this and the actual number of changepoints
        comparison_file.write(str(actual_k) + "\n")
        dim_mean_val = 0
        for i in range(len(full_k_dist)):
            dim_mean_val += float(i) * full_k_dist[i]
        comparison_file.write(str(dim_mean_val) + "\n")

        comparison_file.write(str(actual_effective_k) + "\n")
        eff_dim_mean_val = 0
        for i in range(len(full_effective_k_dist)):
            eff_dim_mean_val += float(i) * full_effective_k_dist[i]
        comparison_file.write(str(eff_dim_mean_val) + "\n")

        comparison_file.write(str(actual_min_proportion_similarity) + "\n")
        min_prop_mean_val = 0
        for i in range(len(min_proportion_trace)):
            min_prop_mean_val += min_proportion_trace[i]
        min_prop_mean_val /= float(len(min_proportion_trace))
        comparison_file.write(str(min_prop_mean_val) + "\n")
        
        # also write the 5%, 25%, 50%, 75% and 95% quantiles to file
        sorted_min_proportion_trace = sorted(min_proportion_trace)
        comparison_file.write(str(sorted_min_proportion_trace[int(5 * len(sorted_min_proportion_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_min_proportion_trace[int(25 * len(sorted_min_proportion_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_min_proportion_trace[int(50 * len(sorted_min_proportion_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_min_proportion_trace[int(75 * len(sorted_min_proportion_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_min_proportion_trace[int(95 * len(sorted_min_proportion_trace) / 100)]) + "\n")

        comparison_file.write(str(actual_similarity) + "\n")
        similarity_mean_val = 0
        for i in range(len(similarity_trace)):
            similarity_mean_val += similarity_trace[i]
        similarity_mean_val /= float(len(similarity_trace))
        comparison_file.write(str(similarity_mean_val) + "\n")
        
        # also write the 5%, 25%, 50%, 75% and 95% quantiles to file
        sorted_similarity_trace = sorted(similarity_trace)
        comparison_file.write(str(sorted_similarity_trace[int(5 * len(sorted_similarity_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_similarity_trace[int(25 * len(sorted_similarity_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_similarity_trace[int(50 * len(sorted_similarity_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_similarity_trace[int(75 * len(sorted_similarity_trace) / 100)]) + "\n")
        comparison_file.write(str(sorted_similarity_trace[int(95 * len(sorted_similarity_trace) / 100)]))
        
    if os.path.isfile(new_filename + "_full_effective_dimension_distribution.txt"):
        os.system("mv " + new_filename + "_full_effective_dimension_distribution.txt ./full_effective_dimension_distributions")
    if os.path.isfile(new_filename + "_full_dimension_distribution.txt"):
        os.system("mv " + new_filename + "_full_dimension_distribution.txt ./full_dimension_distributions")
    if os.path.isfile(new_filename + "_min_proportion_similarity_matrices.txt"):
        os.system("mv " + new_filename + "_min_proportion_similarity_matrices.txt ./min_proportion_similarity_matrices")
    if os.path.isfile(new_filename + "_similarity_matrices.txt"):
        os.system("mv " + new_filename + "_similarity_matrices.txt ./similarity_matrices")
    os.system("rm " + new_filename + "_*")
    os.system("mv comparison_results_" + new_filename + ".txt ./comparison_results")

os.system("mkdir full_effective_dimension_distributions")
os.system("mkdir full_dimension_distributions")
os.system("mkdir min_proportion_similarity_matrices")
os.system("mkdir similarity_matrices")
os.system("mkdir comparison_results")

#pool = Pool(processes = mp.cpu_count())
pool = Pool(processes = 4)
pool.map(run_simulation, range(1000))
#run_simulation(0)
os.system("Rscript simulation_analysis.R")
