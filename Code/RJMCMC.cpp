#include <iostream>
#include <stdio.h>
#include "multiple_processes_regime.h"
#include "basic_RJMCMC.h"
using namespace std;

int main(int argc, char * argv[]) {
    
    string data_file = "";
    string time_file = "uniform_arrivals";
    string model_type = "gaussian";
    unsigned long int start_time = 0;
    unsigned long int intercept = 0;
    unsigned long int end_time = 1;
    unsigned long int diff = 1;
    double poisson_alpha = 1;
    double poisson_beta = 1;
    double markov_alpha = 1;
    double gaussian_mu_0 = 0;
    double gaussian_kappa_0 = 1;
    double gaussian_alpha_0 = 0;
    double gaussian_beta_0 = 0;
    double p = 1;
    double var_p = 1;
    double beta_alpha = 1;
    double dirichlet_alpha = 1;
    double rho = 0.9;
    unsigned long int burnin = 0;
    unsigned long int iterations = 1;
    unsigned long int thinning = 1;
    unsigned long int number_of_association_matrix_bins = 0;
    vector< unsigned long int > starting_changepoints = vector< unsigned long int >(0);
    vector< unsigned long int > separators = vector< unsigned long int >(0);
    vector< unsigned long int > trace_lengths = vector< unsigned long int >(0);
    string starting_changepoints_file = "no_starting_changepoints"; // feed in any non-separator changepoints for the particle to have at the beginning
    string separator_file = "no_separators";
    string trace_lengths_file = "no_trace_lengths";
    unsigned int seed = 2314;
    bool record_basic_rj = false;
    bool record_binary_rj = false;
    bool record_full_rj = true;
    bool record_similarity_matrix = false;
    
    if(argc > 1){
        ifstream InputParameters(argv[1], ios::in);
        InputParameters >> data_file;
        InputParameters >> end_time;
        InputParameters >> diff;
        InputParameters >> gaussian_mu_0;
        InputParameters >> gaussian_kappa_0;
        InputParameters >> gaussian_alpha_0;
        InputParameters >> gaussian_beta_0;
        InputParameters >> poisson_alpha;
        InputParameters >> poisson_beta;
        InputParameters >> markov_alpha;
        InputParameters >> p;
        InputParameters >> var_p;
        InputParameters >> beta_alpha;
        InputParameters >> dirichlet_alpha;
        InputParameters >> rho;
        InputParameters >> burnin;
        InputParameters >> iterations;
        InputParameters >> thinning;
        InputParameters >> number_of_association_matrix_bins;
        InputParameters >> record_basic_rj;
        InputParameters >> record_binary_rj;
        InputParameters >> record_full_rj;
        InputParameters >> record_similarity_matrix;
        InputParameters >> starting_changepoints_file;
        InputParameters >> separator_file;
        InputParameters >> trace_lengths_file;
        InputParameters >> seed;
    }
    
    std::vector< double > model_parameters;
    model_parameters.push_back(gaussian_mu_0);
    model_parameters.push_back(gaussian_kappa_0);
    model_parameters.push_back(gaussian_alpha_0);
    model_parameters.push_back(gaussian_beta_0);
    model_parameters.push_back(poisson_alpha);
    model_parameters.push_back(poisson_beta);
    model_parameters.push_back(markov_alpha);
    
    cout << "data file: " << data_file << endl;
    cout << "end time: " << end_time << endl;
    cout << "diff: " << diff << endl;
    cout << "Gaussian mu_0: " << gaussian_mu_0 << endl;
    cout << "Gaussian kappa_0: " << gaussian_kappa_0 << endl;
    cout << "Gaussian alpha_0: " << gaussian_alpha_0 << endl;
    cout << "Gaussian beta_0: " << gaussian_beta_0 << endl;
    cout << "Poisson alpha: " << poisson_alpha << endl;
    cout << "Poisson beta: " << poisson_beta << endl;
    cout << "Markov alpha: " << markov_alpha << endl;
    cout << "Expected number of change points p: " << p << endl;
    cout << "Variance of p: " << var_p << endl;
    cout << "Binary change points alpha: " << beta_alpha << endl;
    cout << "Regime change points alpha: " << dirichlet_alpha << endl;
    cout << "Rho: " << rho << endl;
    cout << "Burn-in: " << burnin << endl;
    cout << "Iterations: " << iterations << endl;
    cout << "Thinning: " << thinning << endl;
    cout << "Number of association matrix bins: " << number_of_association_matrix_bins << endl;
    cout << "Recording basic change point distns: " << record_basic_rj << endl;
    cout << "Recording binary change point distns: " << record_binary_rj << endl;
    cout << "Recording regime change point distns: " << record_full_rj << endl;
    cout << "Recording similarity matrix: " << record_similarity_matrix << endl;
    cout << "Seed: " << seed << endl;
    
    if (starting_changepoints_file != "no_starting_changepoints") {
        ifstream starting_changepointsFile(starting_changepoints_file);
        unsigned long int new_changepoint;
        while (starting_changepointsFile >> new_changepoint) {
            starting_changepoints.push_back(new_changepoint);
        }
    }
    else {
        cout << "no starting change points" << endl;
    }
    
    if (separator_file != "no_separators") {
        cout << "separators: ";
        ifstream separatorFile(separator_file);
        unsigned long int new_separator;
        while (separatorFile >> new_separator) {
            separators.push_back(new_separator);
            cout << new_separator << '\t';
        }
        cout << endl;
    }
    else {
        cout << "no separators" << endl;
    }
    
    if (trace_lengths_file != "no_trace_lengths") {
        cout << "trace lengths: ";
        ifstream trace_lengthsFile(trace_lengths_file);
        unsigned long int new_trace_length;
        while (trace_lengthsFile >> new_trace_length) {
            trace_lengths.push_back(new_trace_length);
            cout << new_trace_length << '\t';
        }
        cout << endl;
    }
    else {
        cout << "no trace lengths provided, will estimate these from the data" << endl;
    }
    
    mult_process pm = mult_process(end_time, diff, data_file, model_parameters, separators);
    //pm.print_data();
    cout << "finished reading data" << endl;
    
    end_time /= diff;
    p = p / static_cast< double >(end_time);
    for (unsigned int i = 0; i < separators.size(); i++) {
        separators[i] /= diff;
    }
    for (unsigned int i = 0; i < starting_changepoints.size(); i++) {
        starting_changepoints[i] /= diff;
    }
    data_file.erase(data_file.end() - 15, data_file.end());
    
    unsigned long int basic_burnin = burnin;
    unsigned long int basic_iterations = record_basic_rj ? iterations : 0;
    unsigned long int basic_thinning = thinning;
	
    rj rjobject = rj(start_time, end_time, p, var_p, basic_burnin, basic_iterations, basic_thinning, number_of_association_matrix_bins, starting_changepoints, separators, trace_lengths, seed, &pm, intercept, diff);
    rjobject.record_basic_samples(record_basic_rj);
    rjobject.run_basic_simulation();
    if (record_basic_rj) {
        string MAP_cps_Filename = data_file + "_basic_MAP_CPs.txt";
        rjobject.write_basic_MAP_changepoints_to_file(MAP_cps_Filename);
        string dimension_distribution_Filename = data_file + "_basic_dimension_distribution.txt";
        rjobject.write_basic_dimension_distribution_to_file(dimension_distribution_Filename);
        string changepoints_distribution_Filename = data_file + "_basic_changepoints_distribution.txt";
	rjobject.write_basic_changepoints_distribution_to_file(changepoints_distribution_Filename);
   	string log_posterior_trace_Filename = data_file + "_basic_log_posterior_trace.txt";
    	rjobject.write_basic_log_posterior_trace_to_file(log_posterior_trace_Filename);
	string basic_dimension_trace_Filename = data_file + "_basic_dimension_trace.txt";
        rjobject.write_basic_dimension_trace_to_file(basic_dimension_trace_Filename);
    }

    unsigned long int binary_burnin = burnin;
    unsigned long int binary_iterations = record_binary_rj ? iterations : 0;
    unsigned long int binary_thinning = thinning;

    rjobject.set_binary_burnin_iterations_thinning(binary_burnin, binary_iterations, binary_thinning);
    rjobject.convert_basic_particle_to_binary_particle(beta_alpha);
    rjobject.record_binary_samples(record_binary_rj);
    rjobject.run_binary_simulation();
    if (record_binary_rj) {
    	string MAP_cps_Filename = data_file + "_binary_MAP_CPs.txt";
    	rjobject.write_binary_MAP_changepoints_to_file(MAP_cps_Filename);
    	string dimension_distribution_Filename = data_file + "_binary_dimension_distribution.txt";
    	rjobject.write_binary_dimension_distribution_to_file(dimension_distribution_Filename);
	string changepoints_distribution_Filename = data_file + "_binary_changepoints_distribution.txt";
    	rjobject.write_binary_changepoints_distribution_to_file(changepoints_distribution_Filename);
    	string log_posterior_trace_Filename = data_file + "_binary_log_posterior_trace.txt";
    	rjobject.write_binary_log_posterior_trace_to_file(log_posterior_trace_Filename);
    }

    unsigned long int full_burnin = burnin;
    unsigned long int full_iterations = record_full_rj ? iterations : 0;
    unsigned long int full_thinning = thinning;

    rjobject.set_full_burnin_iterations_thinning(full_burnin, full_iterations, full_thinning);
    rjobject.convert_binary_particle_to_full_particle(dirichlet_alpha, rho);
    rjobject.record_full_samples(record_full_rj, data_file);
    rjobject.run_full_simulation();
    if (record_full_rj) {
        string MAP_cps_Filename = data_file + "_full_MAP_CPs.txt";
        rjobject.write_full_MAP_changepoints_to_file(MAP_cps_Filename);
        string dimension_distribution_Filename = data_file + "_full_dimension_distribution.txt";
        rjobject.write_full_dimension_distribution_to_file(dimension_distribution_Filename);
        string effective_dimension_distribution_Filename = data_file + "_full_effective_dimension_distribution.txt";
        rjobject.write_full_effective_dimension_distribution_to_file(effective_dimension_distribution_Filename);
        string changepoints_distribution_Filename = data_file + "_full_changepoints_distribution.txt";
        rjobject.write_full_changepoints_distribution_to_file(changepoints_distribution_Filename, full_iterations);
        string number_of_regimes_Filename = data_file + "_full_number_of_regimes.txt";
        rjobject.write_number_of_regimes_to_file(number_of_regimes_Filename);
        string number_of_observed_regimes_Filename = data_file + "_full_number_of_observed_regimes.txt";
        rjobject.write_number_of_observed_regimes_to_file(number_of_observed_regimes_Filename);
        string log_posterior_trace_Filename = data_file + "_full_log_posterior_trace.txt";
        rjobject.write_full_log_posterior_trace_to_file(log_posterior_trace_Filename);
        string dimension_trace_Filename = data_file + "_full_dimension_trace.txt";
        rjobject.write_dimension_trace_to_file(dimension_trace_Filename);
        string number_of_regimes_trace_Filename = data_file + "_number_of_regimes_trace.txt";
        rjobject.write_number_of_regimes_trace_to_file(number_of_regimes_trace_Filename);
        string full_acceptance_probabilities_Filename = data_file + "_full_acceptance_probabilities.txt";
        rjobject.write_full_acceptance_probabilities_to_file(full_acceptance_probabilities_Filename);
    }
    if (record_similarity_matrix) {
        string similarity_Filename = data_file + "_similarity_matrix.txt";
        rjobject.write_similarity_matrix_to_file(similarity_Filename);
        string min_proportion_similarity_Filename = data_file + "_min_proportion_similarity_matrix.txt";
        rjobject.write_min_proportion_similarity_matrix_to_file(min_proportion_similarity_Filename);
	string similarity_matrices_Filename = data_file + "_similarity_matrices.txt";
	rjobject.write_similarity_matrices_to_file(similarity_matrices_Filename);
	string min_proportion_similarity_matrices_Filename = data_file + "_min_proportion_similarity_matrices.txt";
	rjobject.write_min_proportion_similarity_matrices_to_file(min_proportion_similarity_matrices_Filename);
    }
    if (number_of_association_matrix_bins > 0) { // if recording the association matrix
        string association_matrix_Filename = data_file + "_association_matrix.txt";
        rjobject.write_association_matrix_to_file(association_matrix_Filename);
    }

    return 0;
}
