#ifndef BASIC_PARTICLE_H
#define BASIC_PARTICLE_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "basic_changepoint.h"
#include "binary.h"
#include "regime.h"

class particle{
public:
	particle(const unsigned long int & start_time = 0, const unsigned long int & end_time = 1, const vector< unsigned long int > & separators = vector< unsigned long int >(0), const changepoint & intercept_changepoint = changepoint(), const double & p = 1, const double & var_p = 0); // this sets m_dimension = 0, sets if m_random_p is true, sets m_var_p
	void increase_log_likelihood(const double & likelihood_change);
	void increase_log_k_prior(const double & log_k_prior_change){ m_log_k_prior += log_k_prior_change; }
	void increase_log_binary_I_prior(const double & log_binary_I_prior_change) { m_log_binary_I_prior += log_binary_I_prior_change; }
	void increase_log_full_I_prior(const double & log_full_I_prior_change, const bool & adding_new_regime, const unsigned int & process);
    void increase_log_full_I_prior(const double & log_full_I_prior_change, const bool & adding_new_regime, const unsigned int & process, const unsigned int & trace_index, const bool & adding_new_trace);
	void decrease_log_full_I_prior(const double & log_full_I_prior_change, const bool & removing_unobserved_regime, const bool & adding_new_regime, const unsigned int & process);
	void increase_log_full_I_prior_unobserved(const double & log_full_I_prior_ratio, const vector< int > & unobserved_regime_change);
	void increase_log_regimes_prior(const double & log_regimes_prior) { m_log_regimes_prior += log_regimes_prior; }
	void set_log_likelihood(const double & log_likelihood) { m_log_likelihood = log_likelihood; }
	void set_changepoint_log_likelihood(const int & index, const double & log_likelihood);
	void calculate_and_set_log_k_prior();
	void calculate_and_set_log_binary_I_prior(const size_t & number_of_processes);
	double calculate_and_get_basic_log_likelihood(const size_t & number_of_processes);
	double calculate_and_get_binary_log_likelihood(const size_t & number_of_processes);
	double calculate_and_get_log_binary_I_prior(const size_t & number_of_processes);
	double calculate_and_get_log_full_I_prior(const size_t & number_of_processes);
	void set_log_binary_I_prior(const double & prior) { m_log_binary_I_prior = prior; }
    void set_log_full_I_prior(const double & prior) { m_log_full_I_prior = prior; }
    double calculate_and_get_log_regimes_prior(const size_t & number_of_processes);
    void set_log_regimes_prior(const double & prior) { m_log_regimes_prior = prior; }
    double calculate_and_get_log_full_separators_prior(const size_t & number_of_processes);
    void set_log_full_separators_prior(const double & prior) { m_log_full_separators_prior = prior; }
    void set_arbitrary_extension_density(const double & density) {m_log_arbitrary_extension_density = density;}
    void increase_arbitrary_extension_density(const double & change) {m_log_arbitrary_extension_density += change;}
    double calculate_and_get_full_log_likelihood(const size_t & number_of_processes);
    double dirichlet_ratio(const vector< unsigned int > & transitions);
    double calculate_and_get_basic_log_posterior(const size_t & number_of_processes);
    double calculate_and_get_binary_log_posterior(const size_t & number_of_processes);
    double calculate_and_get_full_log_posterior(const size_t & number_of_processes);
    void print_likelihood() { std::cout << "the likelihood is " << m_log_likelihood << std::endl; }
    void print_sum_of_regime_likelihoods();
    void print_log_k_prior() { cout << m_log_k_prior << endl; }
    double calculate_and_get_log_k_prior();
    void print_log_full_I_prior() { cout << "log_full_I_prior = " << m_log_full_I_prior << endl; }
    void print_calculated_log_full_I_prior();
    void print_log_regimes_prior() { cout << "log_regimes_prior = " << m_log_regimes_prior << endl; }
    void print_calculated_log_regimes_prior();
    void add_separator_changepoint(const changepoint & separator_changepoint, unsigned int index); // this sets m_include_separator = true, sets the separator index, adds the separator cp to tau
    bool get_include_separator() { return m_include_separator; }
    void add_changepoint(const unsigned int & index, const changepoint & new_changepoint);
    void remove_changepoint(const unsigned int & remove_index);
    void move_changepoint(const unsigned int & index, const unsigned long int & position, const double & left_log_likelihood, const double & right_log_likelihood);
    void change_changepoint(const changepoint & changepoint, const unsigned int & change_index);
    void add_binary_changepoint(const unsigned int & index, const changepoint & new_changepoint, const vector< bool > & accept_cp, const vector< double > & left_log_likelihood, const vector< double > & right_log_likelihood);
    void remove_binary_changepoint(const unsigned int & index, const vector< bool > & remove_effective_cp, const vector< double > & merged_log_likelihood);
    void move_binary_changepoint(const unsigned int & index, const unsigned long int & changepoint_position, const vector< bool > & remove_effective_cp, const vector< bool > & accept_cp, const vector< double > & left_log_likelihood, const vector< double > & right_log_likelihood, const vector< double > & merged_log_likelihood);
    void resample_binary_changepoint(const unsigned int & index, const vector< bool > & remove_effective_cp, const vector< bool > & accept_cp, const vector< double > & left_log_likelihood, const vector< double > & right_log_likelihood, const vector< double > & merged_log_likelihood);
	void add_new_binary(const unsigned int & process, const int & cp_index, const double & log_likelihood);
	void add_to_binary(const unsigned int & process, const int & cp_index, const double & combined_log_likelihood);
	unsigned int get_dimension() { return m_dimension; }
    unsigned int calculate_trace_dimension(const unsigned int & trace_index);
    unsigned long int calculate_trace_length(const unsigned int & trace_index);
    unsigned int calculate_and_get_full_effective_dimension();
	bool is_last_changepoint_effective();
    bool is_changepoint_effective_for_process(const unsigned int & cp_index, const unsigned int & process);
    bool does_changepoint_introduce_new_regime_for_process(const unsigned int & cp_idx, const unsigned int & process);
    changepoint & get_changepoint(int cp_index); // if cp_index = -1 returns m_intercept, otherwise m_tau[cp_intercept]
    double get_log_likelihood() { return m_log_likelihood; }
	bool does_changepoint_exist_in_particle(const unsigned long int & changepoint_position, int lower = 0, int upper = -1); // checks if any elements in m_tau has position changepoint_position and if it does not, sets m_add_changepoint_index to the position that the changepoint would be added
    bool is_changepoint_index_separator_index(const unsigned int & changepoint_index);
    unsigned int get_separator_index(const unsigned int & i) { return m_separator_indices[i]; }
    unsigned int get_add_changepoint_index() { return m_add_changepoint_index; }
    double calculate_and_get_add_cp_proposal_ratio(const unsigned long int & trace_length, const unsigned int & trace_index, const bool & basic_changepoint);
	double calculate_and_get_remove_cp_proposal_ratio(const unsigned long int & trace_length = 0, const unsigned int & trace_index = 0, const bool & basic_changepoint = true);
    double calculate_and_get_add_cp_k_prior_ratio();
    double calculate_and_get_remove_cp_k_prior_ratio();
    vector< double > calculate_and_get_binary_log_I_prior_add_ratio(const unsigned int & process);
    vector< double > calculate_and_get_binary_log_I_prior_remove_ratio(const unsigned int & process, const bool & removing_effective_cp);
	vector< double > calculate_and_get_binary_log_I_prior_move_ratio(const unsigned int & j, const bool & removing_effective_cp);
	double full_log_I_prior_add_ratio(const int & index, const bool & new_regime, const unsigned int & process, const unsigned int & previous_regime, const unsigned int & proposed_regime, const unsigned int & subsequent_regime = -1);
	double full_log_I_prior_remove_ratio(const int & index, const unsigned int & process, const unsigned int & previous_regime, const unsigned int & proposed_regime, const bool & new_regime, const unsigned int & actual_regime, const bool & removing_unobserved_regime, const unsigned int & subsequent_regime = -1);
    double log_resampling_separator_changepoint_prior_ratio(const unsigned int & process, const bool & no_following_regime, const bool & new_regime, const bool & removing_unobserved_regime, const unsigned int & actual_regime = -1, const unsigned int & following_regime = -1, const unsigned int & proposed_regime = -1);
    double calculate_and_get_add_unobserved_regimes_full_I_prior_ratio(const vector< bool > & accepted, const size_t & number_of_processes);
    double calculate_and_get_add_unobserved_regimes_full_I_prior_ratio(const unsigned int & process);
    double calculate_and_get_add_unobserved_regimes_proposal_ratio(const vector< bool > & accepted, const size_t & number_of_processes);
    double calculate_and_get_add_unobserved_regimes_regimes_prior_ratio(const vector< bool > & m_adding_unobserved_regimes, const size_t & number_of_processes);
    double calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio(const vector< bool > & rejected, const size_t & number_of_processes);
    double calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio(const unsigned int & process);
    double calculate_and_get_remove_unobserved_regimes_proposal_ratio(const vector< bool > & rejected, const size_t & number_of_processes);
	double calculate_and_get_remove_unobserved_regimes_regimes_prior_ratio(const vector< bool > & removing_unobserved_regimes, const size_t & number_of_processes);
	vector< double > calculate_and_get_adding_binary_log_I_prior_ratio(const unsigned int & process, const bool & new_trace, const unsigned int & trace_index, const int & new_left_cp_index);
	double calculate_and_get_full_adding_binary_log_I_prior_ratio(const unsigned int & process, const unsigned int & regime, const int & number_of_changepoints, const bool & new_trace, const unsigned int & trace_index, const unsigned int & previous_regime);
    void set_binary_marked_vectors(const size_t & number_of_processes, const vector< vector< double > > & log_likelihoods);
	void initiate_binaries(const size_t & number_of_processes) { m_binaries = vector< vector< binary > >(number_of_processes, vector< binary >(0)); }
	void set_all_binary_marked_vectors_equal_to_0_vectors(const size_t & number_of_processes);
	void add_end_binaries(const size_t & number_of_processes);
    void set_binary_indices_from_full_indices(particle & full_particle, const size_t & number_of_processes);
    void set_binary_log_likelihoods(const vector< vector< double > > & log_likelihoods, const size_t & number_of_processes);
    void add_full_changepoint(const unsigned int & index, const changepoint & new_changepoint, const vector< unsigned int > & new_regimes, const vector< vector< double > > & sufficient_statistics, const vector< vector< double > > & log_likelihoods, const vector< double > & previous_log_likelihoods, const vector< double > & number_of_observations);
    void remove_full_changepoint(const unsigned int & index, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_rev, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_rev, const vector< double > & number_of_observations_rev, const vector< bool > & removing_unobserved_regimes);
	void move_full_changepoint(const unsigned int & index, const changepoint & new_changepoint, const vector< unsigned int > & new_regimes, const bool & tau_h_greater_than_tau_h_prime, const vector< vector< double > > & right_sufficient_statistics, const vector< double > & right_number_of_observations, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & middle_sufficient_statistics, const vector< double > & middle_number_of_observations, const vector< double > & previous_log_likelihoods_without_right_sufficient_statistics, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & previous_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & previous_log_likelihoods_without_middle_sufficient_statistics, const vector< vector< double > > & log_likelihoods_with_right_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes);
    void resample_full_changepoint(const unsigned int & index, const vector< unsigned int > & new_regimes, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & regime_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes);
    void alter_unobserved_regimes(const vector< int > & altering_unobserved_regimes, const size_t & number_of_processes);
	vector< unsigned long int > calculate_and_get_changepoint_histogram(const unsigned long int & m_number_of_changepoint_bins);
    // for binary marked vectors, we record a set of bins for each process. Record how many times each interval contains at least 1 effective changepoint for that process.
    vector< unsigned long int > calculate_and_get_binary_changepoint_histogram(const unsigned long int & number_of_changepoint_bins, const size_t & number_of_processes);
	// for regimes, record a set of bins for each process. Record how many times each interval contains at least 1 effective changepoint for that process.
	vector< unsigned long int > calculate_and_get_full_changepoint_histogram(const unsigned long int & number_of_changepoint_bins, const size_t & number_of_processes);
    void calculate_and_add_similarity_matrices(vector< vector< double > > & similarity_matrix, const size_t & number_of_processes);
    void calculate_and_add_min_proportion_similarity_matrices(vector< vector< double > > & min_proportion_similarity_matrix, const size_t & number_of_processes, const vector< unsigned long int > & actual_number_of_observations);
    void add_to_association_matrix(vector< vector< double > > & association_matrix, const unsigned int & process);
	double get_basic_log_posterior();
    double get_binary_log_posterior();
    double get_full_log_posterior();
    double get_full_log_SMC_posterior();
    double get_log_full_I_prior() { return m_log_full_I_prior; }
    double get_log_regimes_prior() { return m_log_regimes_prior; }
    vector< unsigned long int > get_vector_of_binary_left_positions(const unsigned int process);
    vector< int > get_vector_of_binary_left_indices(const unsigned int process);
    int get_binary_left_index(const size_t & process, const unsigned int & index);
    int get_binary_right_index(const size_t & process, const unsigned int & index);
    unsigned int get_binary_I_i_j(const unsigned int & i, const unsigned int & j);
    unsigned int get_full_I_i_j(const unsigned int & i, const unsigned int & j);
    vector< int > get_right_changepoint_indices(const unsigned int & process, const unsigned int & regime) { return m_regimes[process][regime].get_right_changepoint_indices(); }
    unsigned int get_number_of_right_transitions(const unsigned int & process, const unsigned int & regime) { return m_regimes[process][regime].get_number_of_right_transitions(); }
    double get_binary_log_likelihood(const unsigned int & j, const unsigned int & index);
    void set_beta_alpha(const double & alpha) { m_beta_alpha = alpha; }
    void set_full_marked_vectors(const size_t & number_of_processes, const vector< vector< vector< double > > > & sufficient_statistics, const vector< vector< double > > & number_of_observations);
    void set_full_marked_vectors_without_binaries(const size_t & number_of_processes, const vector< vector< vector< double > > > & sufficient_statistics, const vector< vector< double > > & number_of_observations, const vector< vector< double > > & log_likelihoods);
    void set_all_full_marked_vectors_equal_to_binary_marked_vectors(const size_t & number_of_processes);
    void set_all_regimes_to_be_observed(const size_t & number_of_processes);
    void set_number_of_unobserved_regimes(const vector< unsigned int > & number_of_unobserved_regimes) { m_number_of_unobserved_regimes = number_of_unobserved_regimes; }
    void initiate_regime_vectors(const size_t & number_of_processes) { m_regimes = vector< vector< regime > >(number_of_processes, vector< regime >(0)); }
    void add_new_regime(const unsigned int & process, vector< int > & right_changepoint_indices, vector< unsigned int > & right_transitions, vector< unsigned int > & right_transitions_histogram, vector< double > & sufficient_statistics, double & number_of_observations, double & log_likelihood, size_t & number_of_traces, const bool & new_trace, const unsigned int & trace_index, const unsigned int & previous_regime);
    void add_binary_to_regime(const unsigned int & process, const unsigned int & regime_index, const vector< int > & right_changepoint_indices, const vector< unsigned int > & right_transitions, const vector< unsigned int > & right_transitions_histogram, const vector< double > & extra_sufficient_statistics, const unsigned int & extra_number_of_observations, const double & log_likelihood_with_right_sufficient_statistics, const size_t & number_of_traces, const bool & new_trace, const unsigned int & trace_index, const unsigned int & previous_regime);
    void set_dirichlet_alpha(const double & alpha) { m_dirichlet_alpha = alpha; }
    double get_dirichlet_alpha() { return m_dirichlet_alpha; }
    void set_rho(const double & rho) { m_rho = rho; }
    double get_rho() { return m_rho; }
    size_t get_number_of_binaries(const unsigned int process) { return m_binaries[process].size(); }
	size_t get_number_of_regimes(const unsigned int process) { return m_regimes[process].size(); }
    vector< unsigned int > get_number_of_unobserved_regimes() { return m_number_of_unobserved_regimes; }
    bool is_regime_unobserved(const unsigned int & process, const size_t & regime) { return m_unobserved_regimes[process][regime]; }
    unsigned int get_number_of_observed_regimes(const unsigned int & process);
    bool removing_full_changepoint_leaves_highest_regime_unobserved(const unsigned int & process, const unsigned int & regime);
	bool removing_full_changepoint_leaves_regime_unobserved(const unsigned int & process, const unsigned int & regime);
    vector< bool > any_unobserved_regimes();
    unsigned int get_previous_regime(const int & cp_index, const unsigned int & process);
	vector< double > get_sufficient_statistics(const unsigned int & process, const unsigned int & regime_index);
	double get_regime_log_likelihood(const unsigned int & process, const unsigned int & regime_index);
    bool deleting_left_changepoint_removes_regime(const unsigned int & process, const unsigned int & regime_index);
    void increase_separator_indices_greater_or_equal_to_index(const unsigned int & index);
    void decrease_separator_indices_greater_than_index(const unsigned int & index);
    void check_unobserved_regimes(const unsigned int process);
    void check_separator_changepoints();
    void check_observations_in_traces(const unsigned long int & time);
    void check_regimes_and_right_changepoints();
    void check_full_log_posterior(const bool & always_new_regime);
    void check_transitions_out();
    // smc-only functions
    double calculate_and_get_add_cp_proposal_ratio(const unsigned long int & available_cp_positions, const unsigned long int & total_cp_positions);
    double calculate_and_get_remove_cp_proposal_ratio(const unsigned long int & available_cp_positions, const unsigned long int & total_cp_positions);
    void extend_regime(const unsigned int & process, const unsigned int & previous_regime, const vector< double > & previous_regime_new_sufficient_statistics, const double & previous_regime_new_log_likelihood, const unsigned int & trace_index, const double & number_of_observations);
	void add_full_separator_changepoint(changepoint & separator_changepoint, const vector< unsigned int > & new_regimes, const vector< vector< double > > & sufficient_statistics, const vector< vector< double > > & log_likelihoods, const vector< double > & number_of_observations);
	void add_guaranteed_full_changepoint(changepoint & adding_changepoint, const vector< unsigned int > & new_regimes, const vector< vector< double > > & right_sufficient_statistics, const vector< vector< double > > & log_likelihoods_with_right_sufficient_statistics, const vector< double > & right_number_of_observations);
	void increase_log_separator_prior(const double & full_log_prior, const bool & adding_new_regime, const unsigned int & process);
	void resample_number_of_unobserved_regimes(vector< double > uniform_rvs, vector< unsigned int > new_numbers_of_unobserved_regimes); // can either pass uniform RVs to choose the number of regimes for each process or pass in the number of unobserved regimes for each process
	void permute_process_regimes(const unsigned int & process, const vector< unsigned int > & permutation_vector);
	bool is_regime_observed_before(const unsigned int & cp_index, const unsigned int & process);
	int find_changepoint_with_same_regime(const unsigned int & cp_index, const unsigned int & process);
    void add_to_association_matrix_up_to_time(vector< vector< double > > & association_matrix, const unsigned int & process, const unsigned long int & time);
    // returns the log of the full I prior ratio if we were adding the changepoint that we are proposing removing with regime proposed_regime. Index gives the position of the changepoint that we are proposing to remove. new_regime gives whether the proposed_regime is a new regime, either in the sense that if the changepoint were deleted it would remove regime proposed_regime, or in the sense that if we were proposing a changepoint here we would be proposing a new regime for it. actual_regime gives the actual regime of the changepoint that we are removing
    double full_log_I_prior_smc_remove_ratio(const int & index, const unsigned int & process, const unsigned int & previous_regime, const unsigned int & proposed_regime, const bool & new_regime, const unsigned int & actual_regime, const bool & removing_unobserved_regime, const unsigned int & subsequent_regime = -1);
    double log_smc_resampling_separator_changepoint_prior_ratio(const unsigned int & process, const bool & no_following_regime, const bool & new_regime, const bool & removing_unobserved_regime, const unsigned int & actual_regime = -1, const unsigned int & following_regime = -1, const unsigned int & proposed_regime = -1);
	void remove_full_changepoint_smc(const unsigned int & index, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_rev, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_rev, const vector< double > & number_of_observations_rev, const vector< bool > & removing_unobserved_regimes);
	void move_full_changepoint_smc(const unsigned int & index, const changepoint & new_changepoint, const vector< unsigned int > & new_regimes, const bool & tau_h_greater_than_tau_h_prime, const vector< vector< double > > & right_sufficient_statistics, const vector< double > & right_number_of_observations, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & middle_sufficient_statistics, const vector< double > & middle_number_of_observations, const vector< double > & previous_log_likelihoods_without_right_sufficient_statistics, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & previous_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & previous_log_likelihoods_without_middle_sufficient_statistics, const vector< vector< double > > & log_likelihoods_with_right_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes);
	void resample_full_changepoint_smc(const unsigned int & index, const vector< unsigned int > & new_regimes, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & regime_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes);
    double calculate_and_get_log_resampling_last_changepoint_regimes_probability();
    double full_log_I_prior_add_ratio_always_new_regime(const unsigned int & process, const bool & new_regime);
    double full_log_I_prior_remove_ratio_always_new_regime(const unsigned int & process, const bool & new_regime, const bool & removing_unobserved_regime);
    double calculate_and_get_add_unobserved_regimes_full_I_prior_ratio_always_new_regime(const unsigned int & process);
    double calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio_always_new_regime(const unsigned int & process);

protected:
    unsigned long int m_start;
    unsigned long int m_end;
    vector< unsigned long int > m_separators;
    size_t m_number_of_traces;
    unsigned int m_trace_index;
    unsigned int m_h_trace_index;
    unsigned long int m_intercept;
    bool m_include_separator = false;
    vector< unsigned int > m_separator_indices;
    double m_log_likelihood;
    double m_log_k_prior;
    double m_log_binary_I_prior;
    double m_log_full_I_prior;
    double m_log_regimes_prior;
    double m_log_full_separators_prior;
    double m_log_arbitrary_extension_density;
    unsigned int m_dimension;
    changepoint m_intercept_changepoint;
    vector< bool > m_intercept_binary_marked_vector;
    changepoint m_separator_changepoint;
    vector< changepoint > m_tau;
    vector< vector< binary > > m_binaries;
    vector< vector< regime > > m_regimes;
    vector< vector< bool > > m_unobserved_regimes;
    vector< unsigned int > m_number_of_unobserved_regimes;
    double m_p;
    double m_var_p;
    double m_p_alpha;
    double m_p_beta;
    bool m_random_p;
    double m_beta_alpha;
    double m_dirichlet_alpha;
    double m_rho;
    unsigned int m_add_changepoint_index;
    const gsl_rng_type * r_type;
    gsl_rng * r;
};

particle::particle(const unsigned long int & start_time, const unsigned long int & end_time, const vector< unsigned long int > & separators, const changepoint & intercept_changepoint, const double & p, const double & var_p):m_start(start_time), m_end(end_time), m_separators(separators), m_number_of_traces(m_separators.size() + 1), m_intercept_changepoint(intercept_changepoint), m_p(p), m_var_p(var_p) {
    m_dimension = 0;
    m_trace_index = 0;
    if (m_var_p > 0) {
        m_random_p = true;
        if (m_var_p > m_p * (1 - m_p)) {
            cerr << "var > mean (1 - mean)";
        }
        m_p_alpha = m_p * (m_p * (1 - m_p) / m_var_p - 1);
        m_p_beta = m_p_alpha * (1 - m_p) / m_p;
    } else {
        m_random_p = false;
    }
    gsl_rng_env_setup();
    r_type = gsl_rng_default;
    r = gsl_rng_alloc(r_type);
    gsl_rng_set(r, 1);
}

void particle::increase_log_likelihood(const double & likelihood_change){
    m_log_likelihood += likelihood_change;
}

void particle::increase_log_full_I_prior(const double & log_full_I_prior_change, const bool & adding_new_regime, const unsigned int & process) {
    m_log_full_I_prior += log_full_I_prior_change;
    if (adding_new_regime) {
        m_log_full_I_prior -= log(1 - m_rho);
        m_log_regimes_prior += log(1 - m_rho);
        double number_of_regimes = static_cast< double >(m_regimes[process].size());
        m_log_full_I_prior -= static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
        m_log_full_separators_prior += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
    }
}

void particle::increase_log_full_I_prior(const double & log_full_I_prior_change, const bool & adding_new_regime, const unsigned int & process, const unsigned int & trace_index, const bool & adding_new_trace) {
    m_log_full_I_prior += log_full_I_prior_change;
    if (adding_new_regime) {
        m_log_full_I_prior -= log(1 - m_rho);
        m_log_regimes_prior += log(1 - m_rho);
    }
    double number_of_regimes = static_cast< double >(m_regimes[process].size());
    if (!adding_new_trace) {
        if (adding_new_regime) {
            m_log_full_separators_prior += static_cast< double >(trace_index) * log(static_cast< double >(number_of_regimes) / static_cast< double >(number_of_regimes + 1));
            m_log_full_I_prior -= static_cast< double >(trace_index) * log(static_cast< double >(number_of_regimes) / static_cast< double >(number_of_regimes + 1));
        }
    }
    else {
        if (adding_new_regime) {
            m_log_full_I_prior -= static_cast< double >(trace_index) * log(1 / static_cast< double >(number_of_regimes + 1));
            m_log_full_I_prior -= static_cast< double >(trace_index - 1) * log(static_cast< double >(number_of_regimes));
            m_log_full_separators_prior += static_cast< double >(trace_index) * log(1 / static_cast< double >(number_of_regimes + 1));
            m_log_full_separators_prior += static_cast< double >(trace_index - 1) * log(static_cast< double >(number_of_regimes));
        }
        else {
            m_log_full_I_prior -= log(1 / static_cast< double >(number_of_regimes));
            m_log_full_separators_prior += log(1 / static_cast< double >(number_of_regimes));
        }
    }
}

void particle::decrease_log_full_I_prior(const double & log_full_I_prior_change, const bool & removing_unobserved_regime, const bool & adding_new_regime, const unsigned int & process) {
    m_log_full_I_prior -= log_full_I_prior_change;
    if (removing_unobserved_regime != adding_new_regime) {
        if (removing_unobserved_regime) {
            m_log_full_I_prior += log(1 - m_rho);
            m_log_regimes_prior -= log(1 - m_rho);
            double number_of_regimes = static_cast< double >(m_regimes[process].size());
            m_log_full_I_prior += static_cast< double >(m_separator_indices.size()) * log((number_of_regimes - 1) / (number_of_regimes));
            m_log_full_separators_prior -= static_cast< double >(m_separator_indices.size()) * log((number_of_regimes - 1) / (number_of_regimes));
        }
        else {
            m_log_full_I_prior -= log(1 - m_rho);
            m_log_regimes_prior += log(1 - m_rho);
            double number_of_regimes = static_cast< double >(m_regimes[process].size());
            m_log_full_I_prior -= static_cast< double >(m_separator_indices.size()) * log((number_of_regimes) / (number_of_regimes + 1));
            m_log_full_separators_prior += static_cast< double >(m_separator_indices.size()) * log((number_of_regimes) / (number_of_regimes + 1));
        }
    }
}

void particle::increase_log_full_I_prior_unobserved(const double & log_full_I_prior_ratio, const vector< int > & unobserved_regime_change) {
    m_log_full_I_prior += log_full_I_prior_ratio;
    for (unsigned int process = 0; process < unobserved_regime_change.size(); process++) {
        if (unobserved_regime_change[process] == 1) {
            double number_of_regimes = static_cast< double >(m_regimes[process].size());
            double s = static_cast< double >(m_separator_indices.size());
            m_log_full_I_prior -= s * log(number_of_regimes / (number_of_regimes + 1));
            m_log_full_separators_prior += s * log(number_of_regimes / (number_of_regimes + 1));
        }
        else if (unobserved_regime_change[process] == -1) {
            double number_of_regimes = static_cast< double >(m_regimes[process].size());
            double s = static_cast< double >(m_separator_indices.size());
            m_log_full_I_prior -= s * log(number_of_regimes / (number_of_regimes - 1));
            m_log_full_separators_prior += s * log(number_of_regimes / (number_of_regimes - 1));
        }
    }
}

void particle::set_changepoint_log_likelihood(const int & index, const double & log_likelihood) {
    if (index < 0) {
        m_intercept_changepoint.set_log_likelihood(log_likelihood);
    }
    else {
        m_tau[index].set_log_likelihood(log_likelihood);
    }
}

void particle::calculate_and_set_log_k_prior(){
    if(m_random_p){
        m_log_k_prior = gsl_sf_lnbeta(static_cast< double >(m_dimension) + m_p_alpha, m_end - static_cast< double >(m_dimension)) - gsl_sf_lnbeta(m_p_alpha, m_p_beta);
    } else {
        m_log_k_prior = static_cast< double >(m_dimension) * log(m_p) + static_cast< double >(m_end - m_dimension) * log(1 - m_p);
    }
}

void particle::calculate_and_set_log_binary_I_prior(const size_t & number_of_processes){
    m_log_binary_I_prior = 0;
    double k = static_cast< double >(m_dimension);
    double seps = static_cast< double >(m_number_of_traces - 1);
    for (unsigned int j = 0; j < number_of_processes; j++){
        double k_j = static_cast< double >(m_binaries[j].size()) - 2;
        m_log_binary_I_prior += gsl_sf_lnbeta(m_beta_alpha + k_j - seps, m_beta_alpha + k - k_j) - gsl_sf_lnbeta(m_beta_alpha, m_beta_alpha); // don't need seps in m_beta_alpha + k - k_j as would add and subtract
    }
}

double particle::calculate_and_get_basic_log_likelihood(const size_t & number_of_processes) {
    double basic_log_likelihood = m_intercept_changepoint.get_log_likelihood();
    for (unsigned int cp_idx = 0; cp_idx < m_dimension; cp_idx++) {
        basic_log_likelihood += m_tau[cp_idx].get_log_likelihood();
    }
    return basic_log_likelihood;
}

// calculate log likelihood by proceeding adding the log likelihood values stored in the binaries
double particle::calculate_and_get_binary_log_likelihood(const size_t & number_of_processes) {
    double log_likelihood = 0;
    for (unsigned int j = 0; j < number_of_processes; j++) {
        for (unsigned int i = 0; i < m_binaries[j].size() - 1; i++) {
            log_likelihood += m_binaries[j][i].get_log_likelihood();
        }
    }
    return log_likelihood;
}

double particle::calculate_and_get_log_binary_I_prior(const size_t & number_of_processes) {
    double log_binary_I_prior = 0;
    double k = static_cast< double >(m_dimension);
    double seps = static_cast< double >(m_number_of_traces - 1);
    for (unsigned int j = 0; j < number_of_processes; j++){
        double k_j = static_cast< double >(m_binaries[j].size()) - 2;
        log_binary_I_prior += gsl_sf_lnbeta(m_beta_alpha + k_j - seps, m_beta_alpha + k - k_j) - gsl_sf_lnbeta(m_beta_alpha, m_beta_alpha);
    }
    return log_binary_I_prior;
}

// use m_regimes to calculate the log prior for the regimes
double particle::calculate_and_get_log_full_I_prior(const size_t & number_of_processes) {
    double log_prior = 0;
    for (unsigned int process = 0; process < number_of_processes; process++) {
        for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
            log_prior += dirichlet_ratio(m_regimes[process][regime].get_right_transitions_histogram());
        }
    }
    return log_prior;
}

// use m_regimes to calculate the log regimes prior
double particle::calculate_and_get_log_regimes_prior(const size_t & number_of_processes) {
    double log_prior = static_cast< double >(number_of_processes) * log(m_rho);
    for (unsigned int j = 0; j < number_of_processes; j++) {
        log_prior += static_cast< double >(m_regimes[j].size() - 1) * log(1 - m_rho);
    }
    return log_prior;
}

// use m_regimes to calculate the prior for the separators (since there are no regime transitions to separators)
double particle::calculate_and_get_log_full_separators_prior(const size_t & number_of_processes) {
    double log_separators_prior = 0;
    for (unsigned int j = 0; j < number_of_processes; j++) {
        double r_j = static_cast< double >(m_regimes[j].size());
        double s = static_cast< double >(m_separators.size());
        log_separators_prior += s * log(1 / r_j);
    }
    return log_separators_prior;
}

// calculate the log likelihood for the whole multivariate process using information stored in m_regimes
double particle::calculate_and_get_full_log_likelihood(const size_t & number_of_processes) {
    double log_likelihood = 0;
    for (unsigned int process = 0; process < number_of_processes; process++) {
        for (unsigned int regime_index = 0; regime_index < m_regimes[process].size(); regime_index++) {
            log_likelihood += m_regimes[process][regime_index].get_log_likelihood();
        }
    }
    return log_likelihood;
}

// given a histogram of transitions out of a regime, calculate the contribution of that regime to the full log regime prior
double particle::dirichlet_ratio(const vector< unsigned int > & transitions) {
    double ratio = 0;
    double transitions_out = 0;
    double alpha_sum = 0;
    // for each regime that can be transitioned to, calculate the gamma ratio G(alpha + k_{u,v}) / G(alpha) and add the log of it to ratio
    for (unsigned int regime = 0; regime < transitions.size(); regime++) {
        double trans = static_cast< double >(transitions[regime]);
        transitions_out += trans;
        ratio += gsl_sf_lngamma(trans + m_dirichlet_alpha) - gsl_sf_lngamma(m_dirichlet_alpha);
        alpha_sum += m_dirichlet_alpha;
    }
    // calculate the ratio G(alpha + ... + alpha) / G(k_{u,1} + alpha + ... + k_{u,r} + alpha) and add log of it to ratio
    ratio += gsl_sf_lngamma(alpha_sum) - gsl_sf_lngamma(alpha_sum + transitions_out);
    return ratio;
}

double particle::calculate_and_get_basic_log_posterior(const size_t & number_of_processes) {
    double log_k_prior = 0;
    log_k_prior = calculate_and_get_log_k_prior();
    
    double basic_log_likelihood;
    basic_log_likelihood = calculate_and_get_basic_log_likelihood(number_of_processes);
    
    return basic_log_likelihood + log_k_prior;
}

double particle::calculate_and_get_binary_log_posterior(const size_t & number_of_processes) {
    double log_k_prior = 0;
    log_k_prior = calculate_and_get_log_k_prior();
    
    double binary_log_I_prior;
    binary_log_I_prior = calculate_and_get_log_binary_I_prior(number_of_processes);
    
    double binary_log_likelihood;
    binary_log_likelihood = calculate_and_get_binary_log_likelihood(number_of_processes);
    
    return binary_log_likelihood + binary_log_I_prior + log_k_prior;
}

double particle::calculate_and_get_full_log_posterior(const size_t & number_of_processes) {
    double log_k_prior = 0;
    log_k_prior = calculate_and_get_log_k_prior();
    
    double full_log_I_prior;
    full_log_I_prior = calculate_and_get_log_full_I_prior(number_of_processes);
    
    double full_log_regime_prior;
    full_log_regime_prior = calculate_and_get_log_regimes_prior(number_of_processes);
    
    double full_log_separators_prior;
    full_log_separators_prior = calculate_and_get_log_full_separators_prior(number_of_processes);
    
    double full_log_likelihood;
    full_log_likelihood = calculate_and_get_full_log_likelihood(number_of_processes);
    
    return full_log_likelihood + full_log_regime_prior + full_log_I_prior + log_k_prior + full_log_separators_prior;
}

void particle::print_sum_of_regime_likelihoods() {
    double sum = 0;
    for (unsigned int j = 0; j < m_regimes.size(); j++) {
        for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
            sum += m_regimes[j][i].get_log_likelihood();
        }
    }
    cout << sum << endl;
}

double particle::calculate_and_get_log_k_prior(){
    if(m_random_p){
        return gsl_sf_lnbeta(static_cast< double >(m_dimension) + m_p_alpha, m_end - static_cast< double >(m_dimension)) - gsl_sf_lnbeta(m_p_alpha, m_p_beta);
    } else {
        return static_cast< double >(m_dimension) * log(m_p) + static_cast< double >(m_end - m_dimension) * log(1 - m_p);
    }
}

void particle::print_calculated_log_regimes_prior() {
    // check if the number of unobserved regimes is right
    for (unsigned int j = 0; j < m_unobserved_regimes.size(); j++) {
        unsigned int number_of_unobserved_regimes = 0;
        for (unsigned int i = 0; i < m_unobserved_regimes[j].size(); i++) {
            number_of_unobserved_regimes += m_unobserved_regimes[j][i] ? 1 : 0;
        }
        if (number_of_unobserved_regimes != m_number_of_unobserved_regimes[j]) {
            cerr << "number of unobserved regimes does not match";
        }
    }
    // calculate the log regimes prior
    double prior = 0;
    for (unsigned int j = 0; j < m_unobserved_regimes.size(); j++) {
        prior += log(m_rho) + static_cast< double >(m_regimes[j].size() - 1) * log(1 - m_rho);
    }
    cout << "regimes prior " << prior << endl;
}

// add a separator changepoint (assumed to be at the end of the set of separator changepoints already added)
void particle::add_separator_changepoint(const changepoint & separator_changepoint, unsigned int index) {
    m_separator_indices.push_back(index);
    m_tau.insert(m_tau.begin() + index, separator_changepoint);
    m_include_separator = true;
    m_dimension++;
}

void particle::add_changepoint(const unsigned int & index, const changepoint & new_changepoint) {
    m_tau.insert(m_tau.begin() + index, new_changepoint);
    if (m_include_separator) {
        // increase by 1 the cp_index of any separator that is greater than or equal to index
        increase_separator_indices_greater_or_equal_to_index(index);
    }
    m_dimension++;
}

void particle::remove_changepoint(const unsigned int & remove_index) {
    m_tau.erase(m_tau.begin() + remove_index);
    if (m_include_separator) {
        // increase by 1 the cp_index of any separator that is greater than or equal to index
        decrease_separator_indices_greater_than_index(remove_index);
    }
    m_dimension--;
}

void particle::move_changepoint(const unsigned int & index, const unsigned long int & position, const double & left_log_likelihood, const double & right_log_likelihood) {
    m_tau[index].set_log_likelihood(right_log_likelihood);
    m_tau[index].set_position(position);
    if (index == 0) {
        m_intercept_changepoint.set_log_likelihood(left_log_likelihood);
    }
    else {
        m_tau[index - 1].set_log_likelihood(left_log_likelihood);
    }
}

void particle::change_changepoint(const changepoint & changepoint, const unsigned int & change_index){
    m_tau[change_index] = changepoint;
}

// index gives the position for the new cp, new_changepoint is a changepoint object with position, accept_cp is a vector of bools giving which processes accept the cp, left_log_likelihood gives the log likelihood for the binary to the left of the new cp, and right_log_likelihood gives the log likelihood that would be created by the new cp if it is accepted for each process
void particle::add_binary_changepoint(const unsigned int & index, const changepoint & new_changepoint, const vector< bool > & accept_cp, const vector< double > & left_log_likelihood, const vector< double > & right_log_likelihood){
    m_tau.insert(m_tau.begin() + index, new_changepoint); //insert new cp in correct place in tau
    m_dimension++; //increase dimension
    size_t number_of_processes = accept_cp.size();
    
    if (m_include_separator) {
        increase_separator_indices_greater_or_equal_to_index(index);
    }
    
    //start by assigning the vector of binary indices to be the same as the binary indices for the previous changepoint.
    if (index == 0) {
        m_tau[index].set_binary_index(m_intercept_changepoint.get_binary_index());
    }
    else {
        m_tau[index].set_binary_index(m_tau[index - 1].get_binary_index());
    }
    for (unsigned int j = 0; j < number_of_processes; j++) {
        unsigned int prev_binary_index;
        // calculate the binary index of the previous changepoint
        if (index == 0) {
            prev_binary_index = 0;
        }
        else {
            prev_binary_index = m_tau[index - 1].get_binary_index_row(j);
        }
        if (accept_cp[j]){ //if process j accepts the new cp
            // set the log likelihood for the previous binary
            m_binaries[j][prev_binary_index].set_log_likelihood(left_log_likelihood[j]);
            // insert a new binary beginning at the new changepoint.
            m_binaries[j].insert(m_binaries[j].begin() + prev_binary_index + 1, binary(right_log_likelihood[j], index));
            // for every cp after (and including) index, increase the binary index for process j
            for (unsigned int cp_index = index; cp_index < m_dimension; cp_index++){
                m_tau[cp_index].increase_binary_index_row(j, 1);
            }
            // for every binary after the one that has just been inserted, increase the left index for process j.
            for (unsigned int binary_index = prev_binary_index + 2; binary_index < m_binaries[j].size(); binary_index++) {
                m_binaries[j][binary_index].increase_left_index(1);
            }
        }
        else { // process j does not affect the new cp
            // for every binary after prev_binary_index, we need to increase the left indent because we have inserted a new changepoint.
            for (unsigned int binary_index = prev_binary_index + 1; binary_index < m_binaries[j].size(); binary_index++) {
                m_binaries[j][binary_index].increase_left_index(1);
            }
        }
    }
}

// index gives the position to remove the changepoint from tau, remove_effective_cp details if the changepoint to be removed affected each process, and merged log likelihood gives the log likelihood for the entire interval for the merged binary if the removed cp did affect the process
void particle::remove_binary_changepoint(const unsigned int & index, const vector< bool > & remove_effective_cp, const vector< double > & merged_log_likelihood) {
    m_tau.erase(m_tau.begin() + index); // erase the chosen cp
    m_dimension--; //reduce the dimension
    size_t number_of_processes = remove_effective_cp.size();
    if (m_include_separator) {
        // decrease by 1 the cp_index of any separator that is greater than or equal to index
        decrease_separator_indices_greater_than_index(index);
    }
    for (unsigned int j = 0; j < number_of_processes; j++){
        unsigned int prev_binary_index;
        // calculate the binary index of the previous changepoint (for process j)
        if (index == 0) {
            prev_binary_index = 0;
        }
        else {
            prev_binary_index = m_tau[index - 1].get_binary_index_row(j);
        }
        if (remove_effective_cp[j]){ // if process j was affected by the removed cp
            // set the log likelihood for the merged binary to be the merged log likelihood
            m_binaries[j][prev_binary_index].set_log_likelihood(merged_log_likelihood[j]);
            m_binaries[j].erase(m_binaries[j].begin() + prev_binary_index + 1);
            // every cp after the deleted one is now in a binary with index one less than before the removal
            for (unsigned int cp_index = index; cp_index < m_dimension; cp_index++){
                m_tau[cp_index].increase_binary_index_row(j, -1);
            }
            // every binary after the previous_binary_index now drops its left index by 1
            for (unsigned int binary_index = prev_binary_index + 1; binary_index < m_binaries[j].size(); binary_index++) {
                m_binaries[j][binary_index].increase_left_index(-1);
            }
        }
        else { // if process j was unaffected by the removed cp
            for (unsigned int binary_index = prev_binary_index + 1; binary_index < m_binaries[j].size(); binary_index++) {
                // need to drop the left index for every binary after prev_binary_index
                m_binaries[j][binary_index].increase_left_index(-1);
            }
        }
    }
}

// index gives the position of the changepoint to move (it will keep this index after it moves), changepoint position gives the position to move the changepoint to, remove_effective_cp states whether the cp that is being moved and resampled affects each process, accept_cp states whether the resampled cp affects each process, left log_likelihood gives the log_likelihood for each process if the new cp affects process j, right_log_likelihood is similar, merged_log_likelihood gives the likelihood for the binary that contains the moved changepoint if the changepoint does not affect process j.
void particle::move_binary_changepoint(const unsigned int & index, const unsigned long int & changepoint_position, const vector< bool > & remove_effective_cp, const vector< bool > & accept_cp, const vector< double > & left_log_likelihood, const vector< double > & right_log_likelihood, const vector< double > & merged_log_likelihood) {
    m_tau[index].set_position(changepoint_position); // change the position of the moved cp
    size_t number_of_processes = remove_effective_cp.size();
    for (unsigned int j = 0; j < number_of_processes; j++){
        unsigned int prev_binary_index;
        // calculate the binary index of the previous cp
        if (index == 0) {
            prev_binary_index = 0;
        }
        else {
            prev_binary_index = m_tau[index - 1].get_binary_index_row(j);
        }
        if (remove_effective_cp[j]) {
            if (accept_cp[j]) { // if the cp that was moved affected process j and the new one also affects process j
                m_binaries[j][prev_binary_index].set_log_likelihood(left_log_likelihood[j]);
                m_binaries[j][prev_binary_index + 1].set_log_likelihood(right_log_likelihood[j]);
            }
            else { // if the cp that was moved affected process j and the new one doesn't affect process j
                m_binaries[j][prev_binary_index].set_log_likelihood(merged_log_likelihood[j]);
                m_binaries[j].erase(m_binaries[j].begin() + prev_binary_index + 1);
                // the old cp was in m_binaies[j][prev_binary_index + 1] but this binary is deleted now, so m_tau[index], ..., m_tau[dimension - 1] should all drop their binary index for process j
                for (unsigned int cp_index = index; cp_index < m_dimension; cp_index++) {
                    m_tau[cp_index].increase_binary_index_row(j, -1);
                }
            }
        }
        else {
            if (accept_cp[j]) { // if the cp that was moved didn't affect process j and the new one does affect process j
                // correct the log likelihood for the prev_binary_index binary object
                m_binaries[j][prev_binary_index].set_log_likelihood(left_log_likelihood[j]);
                // need to create a new binary
                m_binaries[j].insert(m_binaries[j].begin() + prev_binary_index + 1, binary(right_log_likelihood[j], index));
                // m_tau[index] was in binary prev_binary_index. Now m_tau[cp_index] and all indices above it should increase their binary index by 1
                for (unsigned int cp_index = index; cp_index < m_dimension; cp_index++){
                    m_tau[cp_index].increase_binary_index_row(j, 1);
                }
            }
        }
    }
}

// index gives the position of the changepoint to resample, remove_effective_cp states whether the cp that is being resampled affects each process, accept_cp states whether the resampled cp affects each process, left log_likelihood gives the log_likelihood for each process if the resampled cp affects process j, right_log_likelihood is similar, merged_log_likelihood gives the likelihood for the binary that contains the resampled changepoint if it doesn't affect process j.
void particle::resample_binary_changepoint(const unsigned int & index, const vector< bool > & remove_effective_cp, const vector< bool > & accept_cp, const vector< double > & left_log_likelihood, const vector< double > & right_log_likelihood, const vector< double > & merged_log_likelihood) {
    size_t number_of_processes = remove_effective_cp.size();
    for (unsigned int j = 0; j < number_of_processes; j++){
        if (remove_effective_cp[j] != accept_cp[j]) { // if they are equal then there is no change so we don't need to do anything
            unsigned int prev_binary_index;
            // calculate the binary index of the previous cp
            if (index == 0) {
                prev_binary_index = 0;
            }
            else {
                prev_binary_index = m_tau[index - 1].get_binary_index_row(j);
            }
            if (accept_cp[j]) { // accept the resampled cp, didn't accept the cp previously
                m_binaries[j][prev_binary_index].set_log_likelihood(left_log_likelihood[j]);
                // create a new binary starting at the resampled cp
                m_binaries[j].insert(m_binaries[j].begin() + prev_binary_index + 1, binary(right_log_likelihood[j], index));
                for (unsigned int cp_index = index; cp_index < m_dimension; cp_index++){
                    // increase the binary index for the resampled cp and all those above it for process j
                    m_tau[cp_index].increase_binary_index_row(j, 1);
                }
            }
            else if (remove_effective_cp[j]) { // accepted the cp before, don't accept it now
                m_binaries[j][prev_binary_index].set_log_likelihood(merged_log_likelihood[j]);
                // remove the binary object that contained the resampled cp
                m_binaries[j].erase(m_binaries[j].begin() + prev_binary_index + 1);
                for (unsigned int cp_index = index; cp_index < m_dimension; cp_index++) {
                    // redo the binary index for m_tau[index], ..., m_tau[m_dimension - 1]
                    m_tau[cp_index].increase_binary_index_row(j, -1);
                }
            }
        }
    }
}

void particle::add_new_binary(const unsigned int & process, const int & cp_index, const double & log_likelihood) {
	m_binaries[process].push_back(binary(log_likelihood, cp_index));
    
    //m_log_likelihood += log_likelihood; don't need this as done before

	//alter the binary index of the changepoint cp_index
    if (-1 < cp_index) {
        m_tau[cp_index].set_binary_index_row(process, static_cast< unsigned int >(m_binaries[process].size() - 1));
    }
}

void particle::add_to_binary(const unsigned int & process, const int & cp_index, const double & combined_log_likelihood) {
    //m_log_likelihood += combined_log_likelihood - m_binaries[process].back().get_log_likelihood();
    
	m_binaries[process].back().set_log_likelihood(combined_log_likelihood);

	//alter the binary index of the changepoint cp_index
    if (-1 < cp_index) {
        m_tau[cp_index].set_binary_index_row(process, static_cast< unsigned int >(m_binaries[process].size() - 1));
    }
}

changepoint & particle::get_changepoint(int cp_index){
    if(cp_index < 0){
        return m_intercept_changepoint;
    } else {
        return m_tau[cp_index];
    }
}

// returns if a cp with changepoint_position exists in the particle. If it doesn't exist, m_add_changepoint_index is set to the add_changepoint_index
bool particle::does_changepoint_exist_in_particle(const unsigned long int & changepoint_position, int lower, int upper){//lower defaults to 0, upper to -1
    if(m_dimension == 0){
        m_add_changepoint_index = 0;
        return false;
    }
	if (upper == -1) {
		upper = m_dimension - 1;
	}
    while (upper >= lower) {
        unsigned int mid = (lower + upper) / 2;
        unsigned long int mid_position = m_tau[ mid ].get_position();
        if(mid_position == changepoint_position)
            return true;
        else if(mid_position < changepoint_position)
            lower = mid + 1;
        else
            upper = mid - 1;
    }
    if (lower >= m_dimension) {
        lower = m_dimension - 1;
    }
    if (m_tau[ lower ].get_position() > changepoint_position) {
        m_add_changepoint_index = lower;
    } else {
        m_add_changepoint_index = lower + 1;
    }
    return false;
}

// searches for changepoint in the separator indices, returning true if it is a separator index. If it is not a separator index, false is returned and m_trace_index is set to be the index of the trace that contains the changepoint index
bool particle::is_changepoint_index_separator_index(const unsigned int & changepoint_index) {
    if (!m_include_separator || changepoint_index == numeric_limits< unsigned int >::max()) {
        m_trace_index = 0;
        return false;
    } else {
        // run a binary search to check if changepoint_index is contained within m_separator_indices
        int low = 0, high = static_cast< unsigned int >(m_separator_indices.size() - 1);
        while (low <= high) {
            int mid = (low + high) / 2;
            if (m_separator_indices[mid] == changepoint_index) {
                m_trace_index = static_cast< unsigned int >(mid + 1);
                return true;
            }
            else if (m_separator_indices[mid] < changepoint_index) {
                low = mid + 1;
            }
            else {
                high = mid - 1;
            }
        }
        if (low >= m_separator_indices.size()) {
            low = static_cast< int >(m_separator_indices.size() - 1);
        }
        if (changepoint_index < m_separator_indices[low]) {
            m_trace_index = static_cast< unsigned int >(low);
        }
        else {
            m_trace_index = static_cast< unsigned int >(low + 1);
        }
        return false;
    }
}

// probability of adding the cp we did is 1/(end - dimension). Probability of the reverse move (removing the cp we are now adding) is 1 / (dimension + 1)
double particle::calculate_and_get_add_cp_proposal_ratio(const unsigned long int & trace_length = 0, const unsigned int & trace_index = 0, const bool & basic_changepoint = true) {
    unsigned int trace_dimension = calculate_trace_dimension(trace_index);
    double proposal_ratio = log(static_cast< double >(trace_length - trace_dimension));
    proposal_ratio -= log(static_cast< double >(trace_dimension) + 1);
    if (trace_dimension == 0) {
        if (basic_changepoint) {
            proposal_ratio -= log(3);
        }
        else {
            proposal_ratio -= log(4);
        }
    }
    return proposal_ratio;
}

// probability of adding the cp we are deleting is 1 / (end - (dimension - 1)). Probability of removing the cp is 1 / (dimension)
double particle::calculate_and_get_remove_cp_proposal_ratio(const unsigned long int & trace_dimension, const unsigned int & trace_index, const bool & basic_changepoint) {
    unsigned long int trace_length = calculate_trace_length(trace_index);
    double proposal_ratio = log(static_cast< double >(trace_dimension));
    proposal_ratio -= log(static_cast< double >(trace_length) - static_cast< double >(trace_dimension - 1));
	if (trace_dimension == 1) {
        if (basic_changepoint) {
            proposal_ratio += log(3);
        }
        else {
            proposal_ratio += log(4);
        }
    }
	return proposal_ratio;
	
	/*if (m_include_separator) {
        return log(static_cast< double >(m_dimension - m_number_of_traces + 1)) - log(m_end - static_cast< double >(m_dimension) + 1); //log(m_dimension - (m_number_of_traces - 1)) as can't choose to remove a separator
    } else {
        return log(static_cast< double >(m_dimension)) - log(m_end - static_cast< double >(m_dimension) + 1);
    }*/
}

double particle::calculate_and_get_add_cp_k_prior_ratio(){
    if (m_random_p) {
        return log(static_cast< double >(m_dimension) + m_p_alpha) - log(m_end - static_cast< double >(m_dimension) + m_p_beta);
    } else {
        return log(m_p) - log(1 - m_p);
    }
}

double particle::calculate_and_get_remove_cp_k_prior_ratio(){
    if(m_random_p){
        return log(m_end - static_cast< double >(m_dimension) + 1 + m_p_beta) - log(static_cast< double >(m_dimension) - 1 + m_p_alpha);
    } else {
        return log(1 - m_p) - log(m_p);
    }
}

vector< double > particle::calculate_and_get_binary_log_I_prior_add_ratio(const unsigned int & process) {
    vector< double > q = vector< double >(2, 0);
    double k = static_cast< double >(m_dimension);
    double k_j = static_cast< double >(m_binaries[process].size()) - 2;
    double seps = static_cast< double >(m_number_of_traces - 1);
    q[0] = log(m_beta_alpha + k - k_j) - log(m_beta_alpha + m_beta_alpha + k - seps); // don't need to do (-m_number_of_trace + 1) for the first one as we would add and subtract this
    q[1] = log(m_beta_alpha + k_j - seps) - log(m_beta_alpha + m_beta_alpha + k - seps);
    return q;
}

vector< double > particle::calculate_and_get_binary_log_I_prior_remove_ratio(const unsigned int & process, const bool & removing_effective_cp){
    vector< double > q = vector< double >(2, 0);
    double k = static_cast< double >(m_dimension);
    double k_j = static_cast< double >(m_binaries[process].size()) - 2;
    double seps = static_cast< double >(m_number_of_traces - 1);
    if (!removing_effective_cp){
        q[0] = log(m_beta_alpha + (k - 1) - k_j) - log(m_beta_alpha + m_beta_alpha + (k - 1) - seps);
        q[1] = log(m_beta_alpha + k_j - seps) - log(m_beta_alpha + m_beta_alpha + (k - 1) - seps);
    } else {
        q[0] = log(m_beta_alpha + (k - 1) - (k_j - 1)) - log(m_beta_alpha + m_beta_alpha + (k - 1) - seps);
        q[1] = log(m_beta_alpha + (k_j - 1) - seps) - log(m_beta_alpha + m_beta_alpha + (k - 1) - seps);
    }
    return q;
}

vector< double > particle::calculate_and_get_binary_log_I_prior_move_ratio(const unsigned int & process, const bool & removing_effective_cp){
	vector< double > q = vector< double >(2, 0);
    double k = static_cast< double >(m_dimension);
    double k_j = static_cast< double >(m_binaries[process].size()) - 2;
    double seps = static_cast< double >(m_number_of_traces - 1);
	if (!removing_effective_cp){
        q[0] = log(m_beta_alpha + (k - 1) - k_j); // = log(m_beta_alpha + k - (k_j + 1)
		q[1] = log(m_beta_alpha + k_j - seps);
	}
	else {
		q[0] = log(m_beta_alpha + k - k_j); // = log(m_beta_alpha + (k - 1) - (k_j - 1))
		q[1] = log(m_beta_alpha + (k_j - 1) - seps);
	}
	return q;
}

// how many change points are there between the two separators? Does not include the separators, can equal 0.
unsigned int particle::calculate_trace_dimension(const unsigned int & trace_index) {
    int lower_changepoint_index = ((trace_index == 0) ? -1 : static_cast< int >(m_separator_indices[trace_index - 1]));
    int right_changepoint_index = ((trace_index == m_separators.size()) ? static_cast< int >(m_dimension) : static_cast< int >(m_separator_indices[trace_index]));
    return static_cast< unsigned int >(right_changepoint_index - lower_changepoint_index - 1);
}

unsigned long int particle::calculate_trace_length(const unsigned int & trace_index) {
    unsigned long int lower_bound = ((trace_index == 0) ? 0 : m_separators[trace_index - 1]);
    unsigned long int upper_bound = ((trace_index == m_number_of_traces - 1) ? m_end + 1 : m_separators[trace_index]);
    return upper_bound - lower_bound - 1;
}

// calculate the number of changepoints where at least one process is affected by the changepoint
unsigned int particle::calculate_and_get_full_effective_dimension() {
    if (m_dimension == 0) {
        return 0;
    }
    unsigned int effective_dimension = 0;
    size_t number_of_processes = m_tau[0].get_full_index().size();
    // create a regime vector for the intercept changepoint
    vector< unsigned int > prev_regime_vector(number_of_processes, 0);
    for (const auto cp:m_tau) {
        auto new_regime_vector = cp.get_full_index();
        bool any_different = false;
        unsigned int proc_idx = 0;
        while (!any_different && proc_idx < number_of_processes) {
            any_different = new_regime_vector[proc_idx] != prev_regime_vector[proc_idx];
            proc_idx++;
        }
        if (any_different) {
            effective_dimension++;
        }
        prev_regime_vector = new_regime_vector;
    }
    return effective_dimension;
}

bool particle::is_last_changepoint_effective() {
	if (m_dimension == 0) {
		return true;
	}
	size_t number_of_processes = m_tau[0].get_full_index().size();
	vector< unsigned int > final_regime_vector = m_tau.back().get_full_index();
	vector< unsigned int > penultimate_regime_vector;
	if (m_dimension == 1) {
		penultimate_regime_vector = vector< unsigned int >(number_of_processes, 0);
    }
	else {
		penultimate_regime_vector = m_tau[m_dimension - 2].get_full_index();
	}
	bool any_different = false;
	unsigned int proc_idx = 0;
	while (!any_different && proc_idx < number_of_processes) {
		any_different = final_regime_vector[proc_idx] != penultimate_regime_vector[proc_idx];
		proc_idx++;
	}
	return any_different;
}

bool particle::is_changepoint_effective_for_process(const unsigned int & cp_index, const unsigned int & process) {
    if (m_dimension <= cp_index) {
        cerr << "change point index is too high" << endl;
    }
    if (m_regimes.size() <= process) {
        cerr << "process number is too high" << endl;
    }
    
    if (cp_index == 0) {
        return m_tau[0].get_full_index_row(process) != 0;
    }
    else {
        return m_tau[cp_index].get_full_index_row(process) != m_tau[cp_index - 1].get_full_index_row(process);
    }
}

bool particle::does_changepoint_introduce_new_regime_for_process(const unsigned int & cp_idx, const unsigned int & process) {
    if (m_dimension <= cp_idx) {
        cerr << "change point index is too high" << endl;
    }
    if (m_regimes.size() <= process) {
        cerr << "process number is too high" << endl;
    }
    
    unsigned int regime = m_tau[cp_idx].get_full_index_row(process);
    if (regime == 0) {
        return false;
    }
    
    unsigned int first_occurence_of_regime = static_cast< unsigned int >(m_regimes[process][regime].get_right_changepoint_indices()[0] - 1); // know this isn't -1 because we already checked if the regime is 0
    return first_occurence_of_regime == cp_idx;
}

// calculates the prior ratio for the marked vectors for adding a proposed_regime at index for process. If new_regime is true then this is a new regime, if the changepoint has no subsequent regime then leave the last argument empty.
double particle::full_log_I_prior_add_ratio(const int & index, const bool & new_regime, const unsigned int & process, const unsigned int & previous_regime, const unsigned int & proposed_regime, const unsigned int & subsequent_regime) {
    double q;
    double number_of_regimes = static_cast< double >(get_number_of_regimes(process));
	if (subsequent_regime == numeric_limits< unsigned int >::max()) {
		if (new_regime) {
            q = log(m_dirichlet_alpha);
            double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions());
            q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime) - gsl_sf_lngamma(m_dirichlet_alpha + 1 + number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
            q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            for (unsigned int regime = 0; regime < previous_regime; regime++) {
                double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            for (unsigned int regime = previous_regime + 1; regime < number_of_regimes; regime++) {
                double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            q += log(1 - m_rho);
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
		}
		else {
            double transitions_from_previous_to_proposed = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(proposed_regime));
            double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions());
            q = log(m_dirichlet_alpha + transitions_from_previous_to_proposed) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
		}
	}
	else {
		if (new_regime) {
            double transitions_from_previous_to_subsequent_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(subsequent_regime));
            q = 2 * log(m_dirichlet_alpha) - log(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - log(m_dirichlet_alpha + transitions_from_previous_to_subsequent_regime - 1);
            for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) + gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            q += log(1 - m_rho);
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
		}
		else {
            if (proposed_regime == previous_regime) {
                double transitions_from_previous_to_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(previous_regime));
                double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions());
                q = log(m_dirichlet_alpha + transitions_from_previous_to_previous_regime) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
            }
            else {
                if (subsequent_regime == proposed_regime) {
                    double transitions_from_subsequent_to_subsequent_regime = static_cast< double >(m_regimes[process][subsequent_regime].get_right_transitions_histogram_element(subsequent_regime));
                    double transitions_out_of_subsequent_regime = static_cast< double >(m_regimes[process][subsequent_regime].get_number_of_right_transitions());
                    q = log(m_dirichlet_alpha + transitions_from_subsequent_to_subsequent_regime) - log(m_dirichlet_alpha * number_of_regimes + transitions_out_of_subsequent_regime);
                }
                else {
                    double transitions_from_previous_to_proposed_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(proposed_regime));
                    double transitions_from_proposed_to_subsequent_regime = static_cast< double >(m_regimes[process][proposed_regime].get_right_transitions_histogram_element(subsequent_regime));
                    double transitions_from_previous_to_subsequent_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(subsequent_regime));
                    double transitions_out_of_proposed_regime = static_cast< double >(m_regimes[process][proposed_regime].get_number_of_right_transitions());
                    q = log(m_dirichlet_alpha + transitions_from_previous_to_proposed_regime) + log(m_dirichlet_alpha + transitions_from_proposed_to_subsequent_regime) - log(m_dirichlet_alpha + transitions_from_previous_to_subsequent_regime - 1) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_proposed_regime);
                }
            }
		}
	}
    return q;
}

// returns the log of the full I prior ratio if we were adding the changepoint that we are proposing removing with regime proposed_regime. Index gives the position of the changepoint that we are proposing to remove. new_regime gives whether the proposed_regime is a new regime, either in the sense that if the changepoint were deleted it would remove regime proposed_regime, or in the sense that if we were proposing a changepoint here we would be proposing a new regime for it. actual_regime gives the actual regime of the changepoint that we are removing
double particle::full_log_I_prior_remove_ratio(const int & index, const unsigned int & process, const unsigned int & previous_regime, const unsigned int & proposed_regime, const bool & new_regime, const unsigned int & actual_regime, const bool & removing_unobserved_regime, const unsigned int & subsequent_regime) {
    double q;
    // check if the number of regimes has been reduced before running this, if so then reduce the number of regimes by 1.
    double number_of_regimes = static_cast< double >(get_number_of_regimes(process)) - (removing_unobserved_regime ? 1.0 : 0.0);
    if (subsequent_regime == numeric_limits< unsigned int >::max()) {
        if (new_regime) {
            q = log(m_dirichlet_alpha);
            double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions() - 1);
            q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime) - gsl_sf_lngamma(m_dirichlet_alpha + 1 + number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
            q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            for (unsigned int regime = 0; regime < previous_regime; regime++) {
                double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            for (unsigned int regime = previous_regime + 1; regime < (get_number_of_regimes(process) - (removing_unobserved_regime ? 1 : 0)); regime++) {
                double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            q += log(1 - m_rho);
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
        }
        else {
            double transitions_from_previous_to_proposed = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(proposed_regime)) - (proposed_regime == actual_regime ? 1 : 0);
            double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions() - 1);
            q = log(m_dirichlet_alpha + transitions_from_previous_to_proposed) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
        }
    }
    else {
        if (new_regime) {
            double transitions_from_previous_to_subsequent_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(subsequent_regime));
            if (actual_regime == subsequent_regime) {
                if (previous_regime == actual_regime) {
                    transitions_from_previous_to_subsequent_regime--;
                }
            }
            else {
                if (previous_regime != actual_regime) {
                    transitions_from_previous_to_subsequent_regime++;
                }
            }
            q = 2 * log(m_dirichlet_alpha) - log(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - log(m_dirichlet_alpha + transitions_from_previous_to_subsequent_regime - 1);
            for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                transitions_out_of_regime -= (actual_regime == regime ? 1 : 0);
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) + gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            q += log(1 - m_rho);
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
        }
        else {
            if (proposed_regime == previous_regime) {
                double transitions_from_previous_to_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(previous_regime));
                double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions());
                if (previous_regime == actual_regime) {
                    transitions_from_previous_to_previous_regime--;
                    transitions_out_of_previous_regime--;
                }
                else {
                    if (previous_regime == subsequent_regime) {
                        transitions_from_previous_to_previous_regime++;
                    }
                }
                q = log(m_dirichlet_alpha + transitions_from_previous_to_previous_regime) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
            }
            else {
                if (subsequent_regime == proposed_regime) {
                    double transitions_from_subsequent_to_subsequent_regime = static_cast< double >(m_regimes[process][subsequent_regime].get_right_transitions_histogram_element(subsequent_regime));
                    double transitions_out_of_subsequent_regime = static_cast< double >(m_regimes[process][subsequent_regime].get_number_of_right_transitions());
                    if (actual_regime == subsequent_regime) {
                        transitions_out_of_subsequent_regime--;
                        transitions_from_subsequent_to_subsequent_regime--;
                    }
                    else {
                        if (previous_regime == subsequent_regime) {
                            transitions_from_subsequent_to_subsequent_regime++;
                        }
                    }
                    q = log(m_dirichlet_alpha + transitions_from_subsequent_to_subsequent_regime) - log(m_dirichlet_alpha * number_of_regimes + transitions_out_of_subsequent_regime);
                }
                else { // know proposed != previous and proposed != subsequent
                    // therefore these three transition counts are definitely different
                    double transitions_from_previous_to_proposed_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(proposed_regime));
                    double transitions_from_proposed_to_subsequent_regime = static_cast< double >(m_regimes[process][proposed_regime].get_right_transitions_histogram_element(subsequent_regime));
                    double transitions_from_previous_to_subsequent_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(subsequent_regime));
                    if (actual_regime == subsequent_regime) {
                        if (previous_regime == actual_regime) {
                            transitions_from_previous_to_subsequent_regime--;
                        }
                        else {
                            if (proposed_regime == actual_regime) {
                                transitions_from_proposed_to_subsequent_regime--;
                            }
                        }
                    }
                    else {
                        if (previous_regime != actual_regime) {
                            if (proposed_regime == actual_regime) {
                                transitions_from_previous_to_proposed_regime--;
                                transitions_from_proposed_to_subsequent_regime--;
                            }
                            transitions_from_previous_to_subsequent_regime++;
                        }
                    }
                    double transitions_out_of_proposed_regime = static_cast< double >(m_regimes[process][proposed_regime].get_number_of_right_transitions());
                    if (proposed_regime == actual_regime) {
                        transitions_out_of_proposed_regime--;
                    }
                    q = log(m_dirichlet_alpha + transitions_from_previous_to_proposed_regime) + log(m_dirichlet_alpha + transitions_from_proposed_to_subsequent_regime) - log(m_dirichlet_alpha + transitions_from_previous_to_subsequent_regime - 1) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_proposed_regime);
                }
            }
        }
    }
    return q;
}

// following_regime and proposed_regime default to -1. Gives the log of the full I prior ratio for adding the separator as a changepoint after it has been removed.
double particle::log_resampling_separator_changepoint_prior_ratio(const unsigned int & process, const bool & no_following_regime, const bool & new_regime, const bool & removing_unobserved_regime, const unsigned int & actual_regime, const unsigned int & following_regime, const unsigned int & proposed_regime) {
    double q = 0;
    double number_of_regimes = static_cast< double >(get_number_of_regimes(process)) - (removing_unobserved_regime ? 1.0 : 0.0);
    if (no_following_regime) {
        if (new_regime) {
            for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
                double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
            q += log(1 - m_rho);
        }
    }
    else {
        if (new_regime) {
            for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
                double transitions_out_of_regime;
                if (regime == actual_regime) { // one transition has been removed from the actual regime.
                    transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions()) - 1;
                }
                else {
                    transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                }
                q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            }
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
            q += log(m_dirichlet_alpha) - log(number_of_regimes * m_dirichlet_alpha + m_dirichlet_alpha);
            q += log(1 - m_rho);
        }
        else {
            if (actual_regime == proposed_regime) {
                double transitions_from_proposed_to_following_regime = static_cast< double >(m_regimes[process][proposed_regime].get_right_transitions_histogram_element(following_regime));
                double transitions_out_of_proposed_regime = static_cast< double >(m_regimes[process][proposed_regime].get_number_of_right_transitions());
                q += log(transitions_from_proposed_to_following_regime - 1 + m_dirichlet_alpha) - log(transitions_out_of_proposed_regime - 1 + (number_of_regimes * m_dirichlet_alpha));
            }
            else {
                double transitions_from_proposed_to_following_regime = static_cast< double >(m_regimes[process][proposed_regime].get_right_transitions_histogram_element(following_regime));
                double transitions_out_of_proposed_regime = static_cast< double >(m_regimes[process][proposed_regime].get_number_of_right_transitions());
                q += log(transitions_from_proposed_to_following_regime + m_dirichlet_alpha) - log(transitions_out_of_proposed_regime + number_of_regimes * m_dirichlet_alpha);
            }
        }
    }
    return q;
}

double particle::calculate_and_get_add_unobserved_regimes_full_I_prior_ratio(const vector< bool > & accepted, const size_t & number_of_processes) {
    double ratio = 0;
    for (unsigned int process = 0; process < number_of_processes; process++) {
        if (accepted[process]) {
            double number_of_regimes = static_cast< double >(m_regimes[process].size());
            for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
                ratio += gsl_sf_lngamma((number_of_regimes + 1) * m_dirichlet_alpha) + gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions())) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma((number_of_regimes + 1) * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions()));
            }
            ratio += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
        }
    }
    return ratio;
}

double particle::calculate_and_get_add_unobserved_regimes_full_I_prior_ratio(const unsigned int & process) {
    double ratio = 0;
    double number_of_regimes = static_cast< double >(m_regimes[process].size());
    for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
        ratio += gsl_sf_lngamma((number_of_regimes + 1) * m_dirichlet_alpha) + gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions())) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma((number_of_regimes + 1) * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions()));
    }
    ratio += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
    return ratio;
}

double particle::calculate_and_get_add_unobserved_regimes_proposal_ratio(const vector< bool > & accepted, const size_t & number_of_processes) {
    double ratio = 0;
    for (unsigned int process = 0; process < number_of_processes; process++) {
        if (accepted[process]) {
            ratio += log(static_cast< double >(m_regimes[process].size()));
            ratio -= log(static_cast< double >(m_number_of_unobserved_regimes[process]) + 1);
        }
    }
    return ratio;
}

double particle::calculate_and_get_add_unobserved_regimes_regimes_prior_ratio(const vector< bool > & adding_unobserved_regimes, const size_t & number_of_processes) {
    double ratio = 0;
    for (unsigned int j = 0; j < number_of_processes; j++) {
        if (adding_unobserved_regimes[j]) {
            ratio += log(1 - m_rho);
        }
    }
    return ratio;
}

double particle::calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio(const vector< bool > & rejected, const size_t & number_of_processes) {
    double ratio = 0;
    for (unsigned int process = 0; process < number_of_processes; process++) {
        if (rejected[process]) {
            double number_of_regimes = static_cast< double >(m_regimes[process].size());
            for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
                ratio += gsl_sf_lngamma((number_of_regimes - 1) * m_dirichlet_alpha) + gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions())) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma((number_of_regimes - 1) * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions()));
            }
            ratio += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes - 1));
        }
    }
    return ratio;
}

double particle::calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio(const unsigned int & process) {
    double ratio = 0;
    double number_of_regimes = static_cast< double >(m_regimes[process].size());
    for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
        ratio += gsl_sf_lngamma((number_of_regimes - 1) * m_dirichlet_alpha) + gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions())) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma((number_of_regimes - 1) * m_dirichlet_alpha + static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions()));
    }
    ratio += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes - 1));
    return ratio;
}

double particle::calculate_and_get_remove_unobserved_regimes_proposal_ratio(const vector< bool > & rejected, const size_t & number_of_processes) {
    double ratio = 0;
    for (unsigned int process = 0; process < number_of_processes; process++) {
        if (rejected[process]) {
            ratio -= log(static_cast< double >(m_regimes[process].size()) - 1);
            ratio += log(static_cast< double >(m_number_of_unobserved_regimes[process]));
        }
    }
    return ratio;
}

double particle::calculate_and_get_remove_unobserved_regimes_regimes_prior_ratio(const vector< bool > & removing_unobserved_regimes, const size_t & number_of_processes) {
	double ratio = 0;
	for (unsigned int j = 0; j < number_of_processes; j++) {
		if (removing_unobserved_regimes[j]) {
			ratio -= log(1 - m_rho);
		}
	}
	return ratio;
}

vector< double > particle::calculate_and_get_adding_binary_log_I_prior_ratio(const unsigned int & process, const bool & new_trace, const unsigned int & trace_index, const int & new_left_cp_index) {
	if (new_left_cp_index == -1) {
		return vector< double >(2, 0);
	}
	double old_k = static_cast< double >(new_left_cp_index);
	double old_seps = static_cast< double >(trace_index - (new_trace ? 1 : 0));
	double old_k_j = static_cast< double >(m_binaries[process].size()) - 1;
	double old_log_binary_I_prior = gsl_sf_lnbeta(m_beta_alpha + old_k_j - old_seps, m_beta_alpha + old_k - old_k_j) - gsl_sf_lnbeta(m_beta_alpha, m_beta_alpha);

	double new_k = old_k + 1;
	double new_seps = trace_index;
	vector< double > log_binary_prior_ratio(2, 0);
	if (new_trace) { // if we are adding a new trace then we must accept this binary so set the ratio for not changing to -infinity so it won't be chosen
		log_binary_prior_ratio[0] = -1e300; 
	}
	else {
		log_binary_prior_ratio[0] = gsl_sf_lnbeta(m_beta_alpha + old_k_j - new_seps, m_beta_alpha + new_k - old_k_j) - gsl_sf_lnbeta(m_beta_alpha, m_beta_alpha) - old_log_binary_I_prior;
	}
	log_binary_prior_ratio[1] = gsl_sf_lnbeta(m_beta_alpha + (old_k_j + 1) - new_seps, m_beta_alpha + new_k - (old_k_j + 1)) - gsl_sf_lnbeta(m_beta_alpha, m_beta_alpha) - old_log_binary_I_prior;
	return log_binary_prior_ratio;
}

double particle::calculate_and_get_full_adding_binary_log_I_prior_ratio(const unsigned int & process, const unsigned int & regime, const int & number_of_changepoints, const bool & new_trace, const unsigned int & trace_index, const unsigned int & previous_regime) {
    unsigned int number_of_regimes = static_cast< unsigned int >(m_regimes[process].size());
    // set the last transition of the previous regime to be regime (unless we are starting a new trace)
    double log_prior;
    // need to add the number_of_changepoints - 1 transitions from regime -> regime to the transitions_histogram_vector for regime
    if (!new_trace) { // no new trace, so need to add a transition at the end of the previous regime
        if (regime == number_of_regimes) { // need to add 1 to the length of each transition_histogram_vector
            log_prior = 0;
            for (unsigned int i = 0; i < number_of_regimes; i++) {
                vector< unsigned int > transitions_histogram = m_regimes[process][i].get_right_transitions_histogram();
                log_prior -= dirichlet_ratio(transitions_histogram);
                transitions_histogram.push_back(0);
                if (i == previous_regime) {
                    transitions_histogram.back()++;
                }
                log_prior += dirichlet_ratio(transitions_histogram);
            }
            // and now make the transitions_histogram_vector for the new regime
            vector< unsigned int > regime_transitions_histogram(number_of_regimes + 1);
            for (int i = 0; i < number_of_changepoints - 1; i++) {
                regime_transitions_histogram[regime]++;
            }
            log_prior += dirichlet_ratio(regime_transitions_histogram);
            log_prior += log(1 - m_rho);
            if (0 < trace_index) {
                log_prior += static_cast< double >(trace_index) * log(static_cast< double >(number_of_regimes) / static_cast< double >(number_of_regimes + 1));
            }
        }
        else { // log prior is dir_ratio(regime_hist_with_new_transitions) + dir_ratio(prev_hist_with_new_transitions) - dir_ratio(previous_hist) - dir_ratio(regime_hist)
            if (regime != previous_regime) {
                // calculate the log prior ratio for adding these transitions to the regime
                double nr = static_cast< double >(number_of_regimes), kr = static_cast< double >(m_regimes[process][regime].get_transitions_out()), krr = static_cast< double >(m_regimes[process][regime].get_right_transitions_histogram_element(regime));
                log_prior = gsl_sf_lngamma(nr * m_dirichlet_alpha + kr);
                log_prior -= gsl_sf_lngamma(nr * m_dirichlet_alpha + kr + number_of_changepoints - 1);
                log_prior += gsl_sf_lngamma(m_dirichlet_alpha + krr + number_of_changepoints - 1);
                log_prior -= gsl_sf_lngamma(m_dirichlet_alpha + krr);
                /*vector< unsigned int > regime_transitions_histogram = m_regimes[process][regime].get_right_transitions_histogram();
                double flog_prior = -dirichlet_ratio(regime_transitions_histogram);
                // add self-transitions for the binary that we are adding
                for (int i = 0; i < number_of_changepoints - 1; i++) {
                    regime_transitions_histogram[regime]++;
                }
                flog_prior += dirichlet_ratio(regime_transitions_histogram);
                cout << "log prior diff: " << log_prior - flog_prior << endl;*/
                // calculate the log prior ratio for adding the new transition out of the previous regime
                double kpr = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(regime)), kp = static_cast< double >(m_regimes[process][previous_regime].get_transitions_out());
                log_prior += log(m_dirichlet_alpha + kpr);
                log_prior -= log(nr * m_dirichlet_alpha + kp);
            }
            else {
                vector< unsigned int > regime_transitions_histogram = m_regimes[process][regime].get_right_transitions_histogram();
                log_prior = -dirichlet_ratio(regime_transitions_histogram);
                // add self-transitions for the binary that we are adding
                for (int i = 0; i < number_of_changepoints - 1; i++) {
                    regime_transitions_histogram[regime]++;
                }
                // add a transition from the previous regime since this isn't a new trace
                regime_transitions_histogram[regime]++;
                log_prior += dirichlet_ratio(regime_transitions_histogram);
            }
        }
    }
    else { // don't add a transition at the end of the previous_regime
        if (regime == number_of_regimes) { // need to add 1 to the length of each transition_histogram_vector
            log_prior = 0;
            for (unsigned int i = 0; i < number_of_regimes; i++) {
                vector< unsigned int > transitions_histogram = m_regimes[process][i].get_right_transitions_histogram();
                log_prior -= dirichlet_ratio(transitions_histogram);
                transitions_histogram.push_back(0);
                log_prior += dirichlet_ratio(transitions_histogram);
            }
            // and now make the transitions_histogram_vector for the new regime
            vector< unsigned int > regime_transitions_histogram(number_of_regimes + 1);
            for (int i = 0; i < number_of_changepoints - 1; i++) {
                regime_transitions_histogram[regime]++;
            }
            log_prior += dirichlet_ratio(regime_transitions_histogram);
            log_prior += log(1 - m_rho);
            log_prior += static_cast< double >(trace_index) * log(1 / static_cast< double >(number_of_regimes + 1));
            if (0 < trace_index) {
                log_prior += static_cast< double >(trace_index - 1) * log(static_cast< double >(number_of_regimes));
            }
        }
        else { // log prior is dir_ratio(regime_hist_with_new_transitions) + dir_ratio(prev_hist_with_new_transitions) - dir_ratio(previous_hist) - dir_ratio(regime_hist)
            vector< unsigned int > regime_transitions_histogram = m_regimes[process][regime].get_right_transitions_histogram();
            log_prior = -dirichlet_ratio(regime_transitions_histogram);
            for (int i = 0; i < number_of_changepoints - 1; i++) {
                regime_transitions_histogram[regime]++;
            }
            log_prior += dirichlet_ratio(regime_transitions_histogram);
            if (0 < trace_index) {
                log_prior += log(1 / static_cast< double >(number_of_regimes));
            }
        }
    }
    
    return log_prior;
}

// index gives the index for the new changepoint, new_changepoint gives the changepoint object to add, new_regimes is a vector with the new regime to add for each process, sufficient_statistics gives the sufficient statistics in the interval [tau_h^\prime, \tau_h+1), log_likelihoods[j][0] gives the log likelihood for the previous regime and log_likelihoods[j][1] gives the log likelihood for new_regimes[j]
void particle::add_full_changepoint(const unsigned int & index, const changepoint & new_changepoint, const vector< unsigned int > & new_regimes, const vector< vector< double > > & sufficient_statistics, const vector< vector< double > > & log_likelihoods, const vector< double > & previous_log_likelihoods, const vector< double > & number_of_observations) {
    // insert new_changepoint into m_tau
    m_tau.insert(m_tau.begin() + index, new_changepoint);
    m_dimension++; //increase dimension
    size_t number_of_processes = new_regimes.size();
    // assign the regime indices for the new changepoint.
    m_tau[index].set_full_index(new_regimes);
    if (m_include_separator) {
        // increase by 1 the cp_index of any separator that is greater than or equal to index
        increase_separator_indices_greater_or_equal_to_index(index);
    }
    
    for (unsigned int j = 0; j < number_of_processes; j++) {
        unsigned int prev_regime_index;
        // calculate the regime index of the previous changepoint
        if (index == 0) {
            prev_regime_index = 0;
        }
        else {
            prev_regime_index = m_tau[index - 1].get_full_index_row(j);
        }
        // if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
        bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
        if (adding_new_regime) {
            for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
                // increase the size of the transitions histogram for each of the regimes in this process.
                m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
            }
            // insert new regime into m_regimes[j]. We don't include the right changepoint position, right transition, sufficient statistics or number_of_observations here as it will be added later
            // these will be converted to the correct values later
            m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(sufficient_statistics[j].size(), 0), m_number_of_traces, m_trace_index, 0));
            // these unobserved regime values will be corrected later
            m_unobserved_regimes[j].push_back(true);
            m_number_of_unobserved_regimes[j]++;
        }
        
        // increase the values for the right indices in every regime so that index + 1 -> index + 2, ...
        for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
            m_regimes[j][regime].increase_right_indices_greater_than_or_equal_to(index + 1);
        }
        
        // insert the interval and transition out of the new_regime
        // replace the right transition from previous regime to subsequent regime with a transition to new_regime (unless new_regime is equal to subsequent_regime), adding a transition out of previous regime if the new changepoint is inserted at the end of m_tau.
        if (index == m_dimension - 1 || is_changepoint_index_separator_index(index + 1)) {
            // if index + 1 is a separator then m_trace_index is set to be the trace index of index
            m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
            m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);
            if (index == m_dimension - 1) {
                m_trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
            }
        }
        else {
            unsigned int subs_regime_index = m_tau[index + 1].get_full_index_row(j);
            if (new_regimes[j] != subs_regime_index) {
                m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
            }
            m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j], subs_regime_index);
        }
        // calculate the trace index
        is_changepoint_index_separator_index(index); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index + 1) has just been run
        
        // assign right sufficient statistics and observation numbers to the new regime and take them away from the previous regime (unless both regimes are equal). Calculate the log likelihood for each regime
        if (new_regimes[j] != prev_regime_index) {
            m_regimes[j][new_regimes[j]].add_sufficient_statistics(sufficient_statistics[j]);
            m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods[j][new_regimes[j]]);
            m_regimes[j][new_regimes[j]].add_observations(m_trace_index, number_of_observations[j]);
            m_regimes[j][prev_regime_index].remove_sufficient_statistics(sufficient_statistics[j]);
            m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods[j]);
            m_regimes[j][prev_regime_index].remove_observations(m_trace_index, number_of_observations[j]);
        }
    }
}

void particle::remove_full_changepoint(const unsigned int & index, const vector< vector< double > > & right_sufficient_statistics_rev, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_rev, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_rev, const vector< double > & number_of_observations_rev, const vector< bool > & removing_unobserved_regimes) {
    // get the regime indices for the index'th changepoint
    vector< unsigned int > removed_regimes(0);
    removed_regimes = m_tau[index].get_full_index();
    
    // remove the changepoint from m_tau
    m_tau.erase(m_tau.begin() + index);
    m_dimension--;
    size_t number_of_processes = right_sufficient_statistics_rev.size();
    if (m_include_separator) {
        // if there are separators, decrease the index of any separator greater than index
        decrease_separator_indices_greater_than_index(index);
    }
    
    for (unsigned int j = 0; j < number_of_processes; j++) {
        unsigned int prev_regime_index;
        // calculate the regime index of the previous changepoint
        if (index == 0) {
            prev_regime_index = 0;
        }
        else {
            prev_regime_index = m_tau[index - 1].get_full_index_row(j);
        }
        
        // remove the interval and transition out of the removed_regime
        // replace the right transition from previous regime to removed_regime with a transition to subsequent_regime (unless new_regime is equal to subsequent_regime), removing a transition out of previous regime if the changepoint is removed from the end of m_tau.
        // calculate m_trace_index for the changepoint we are removing
        if (is_changepoint_index_separator_index(index) || index == m_dimension) {
            if (index == m_dimension) {
                m_trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
            }
            m_regimes[j][prev_regime_index].remove_right_transition(index);
            m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
        }
        else {
            unsigned int subs_regime_index = m_tau[index].get_full_index_row(j);
            if (removed_regimes[j] != subs_regime_index) {
                m_regimes[j][prev_regime_index].alter_right_transition(index, subs_regime_index);
            }
            m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
        }
        // calculate the trace index
        is_changepoint_index_separator_index(index - 1); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index) has just been run and index now corresponds to the next changepoint
        
        // remove right sufficient statistics from the regime for the deleted cp and add them to the previous regime (unless both regimes are equal). Calculate the log likelihood for each regime
        if (removed_regimes[j] != prev_regime_index) {
            m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_rev[j]);
            m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_log_likelihoods_without_right_sufficient_statistics_rev[j]);
            m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, number_of_observations_rev[j]);
            m_regimes[j][prev_regime_index].add_sufficient_statistics(right_sufficient_statistics_rev[j]);
            m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_right_sufficient_statistics_rev[j]);
            m_regimes[j][prev_regime_index].add_observations(m_trace_index, number_of_observations_rev[j]);
        }
        
        if (removing_unobserved_regimes[j]) {
            // if we are removing the top regime, then this regime must be deleted
            if (!m_unobserved_regimes[j].back()) {
                cerr << "regime to delete is not unobserved" << endl;
            }
            m_regimes[j].erase(m_regimes[j].end() - 1);
            unsigned int regime_to_remove = static_cast< unsigned int >(m_regimes[j].size());
            // for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
            for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
                m_regimes[j][regime].remove_unobserved_regime(regime_to_remove);
            }
            m_unobserved_regimes[j].erase(m_unobserved_regimes[j].end() - 1);
            m_number_of_unobserved_regimes[j]--;
        }
        
        // decrease the values for right indices in every regime so that index + 2 -> index + 1, index + 3 -> index + 2, ...
        for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
            m_regimes[j][regime].decrease_right_indices_greater_than(index + 1);
        }
    }
}

void particle::move_full_changepoint(const unsigned int & index, const changepoint & new_changepoint, const vector< unsigned int > & new_regimes, const bool & tau_h_greater_than_tau_h_prime, const vector< vector< double > > & right_sufficient_statistics, const vector< double > & right_number_of_observations, const vector< vector< double > > & right_sufficient_statistics_reverse,const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & middle_sufficient_statistics, const vector< double > & middle_number_of_observations, const vector< double > & previous_log_likelihoods_without_right_sufficient_statistics, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & previous_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & previous_log_likelihoods_without_middle_sufficient_statistics, const vector< vector< double > > & log_likelihoods_with_right_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes) {
	vector< unsigned int > removed_regimes(0);
	removed_regimes = m_tau[index].get_full_index();

	m_tau[index] = new_changepoint;
	m_tau[index].set_full_index(new_regimes);
	size_t number_of_processes = removed_regimes.size();
	for (unsigned int j = 0; j < number_of_processes; j++) {
		unsigned int prev_regime_index;
		// calculate the regime index of the previous changepoint
		if (index == 0) {
			prev_regime_index = 0;
		}
		else {
			prev_regime_index = m_tau[index - 1].get_full_index_row(j);
		}
        
        // if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
        bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
        if (adding_new_regime) {
            for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
                // increase the size of the transitions histogram for each of the regimes in this process.
                m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
            }
            m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(right_sufficient_statistics[j].size(), 0), m_number_of_traces, m_trace_index, 0));
            // the interval will be added later
            m_unobserved_regimes[j].push_back(true);
            m_number_of_unobserved_regimes[j]++;
        }

		// remove right sufficient statistics from the regime for the deleted cp and add them to the previous regime (unless both regimes are equal). Calculate the log likelihood for each regime
		if (removed_regimes[j] == prev_regime_index) {
			if (prev_regime_index != new_regimes[j]) {
				m_regimes[j][prev_regime_index].remove_sufficient_statistics(right_sufficient_statistics[j]);
				m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_without_right_sufficient_statistics[j]);
                m_regimes[j][prev_regime_index].remove_observations(m_trace_index, right_number_of_observations[j]);
				m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics[j]);
				m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]]);
                m_regimes[j][new_regimes[j]].add_observations(m_trace_index, right_number_of_observations[j]);
			}
		}
		else {
			if (new_regimes[j] == prev_regime_index) {
				m_regimes[j][prev_regime_index].add_sufficient_statistics(right_sufficient_statistics_reverse[j]);
				m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_right_sufficient_statistics_reverse[j]);
                m_regimes[j][prev_regime_index].add_observations(m_trace_index, right_number_of_observations_reverse[j]);
				m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
				m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
                m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, right_number_of_observations_reverse[j]);
			}
			else {
				if (new_regimes[j] == removed_regimes[j]) {
					if (tau_h_greater_than_tau_h_prime) {
						m_regimes[j][prev_regime_index].remove_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_without_middle_sufficient_statistics[j]);
                        m_regimes[j][prev_regime_index].remove_observations(m_trace_index, middle_number_of_observations[j]);
						m_regimes[j][removed_regimes[j]].add_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_with_middle_sufficient_statistics[j]);
                        m_regimes[j][removed_regimes[j]].add_observations(m_trace_index, middle_number_of_observations[j]);
					}
					else {
						m_regimes[j][prev_regime_index].add_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_middle_sufficient_statistics[j]);
                        m_regimes[j][prev_regime_index].add_observations(m_trace_index, middle_number_of_observations[j]);
						m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_middle_sufficient_statistics[j]);
                        m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, middle_number_of_observations[j]);
					}
				}
				else {
					if (tau_h_greater_than_tau_h_prime) {
						m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
                        m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, right_number_of_observations_reverse[j]);
						m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics[j]);
						m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]]);
                        m_regimes[j][new_regimes[j]].add_observations(m_trace_index, right_number_of_observations[j]);
						m_regimes[j][prev_regime_index].remove_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_without_middle_sufficient_statistics[j]);
                        m_regimes[j][prev_regime_index].remove_observations(m_trace_index, middle_number_of_observations[j]);
					}
					else {
						m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
                        m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, right_number_of_observations_reverse[j]);
						m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics[j]);
						m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]]);
                        m_regimes[j][new_regimes[j]].add_observations(m_trace_index, right_number_of_observations[j]);
						m_regimes[j][prev_regime_index].add_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_middle_sufficient_statistics[j]);
                        m_regimes[j][prev_regime_index].add_observations(m_trace_index, middle_number_of_observations[j]);
					}
				}
			}
		}

		// replace the interval for removed_regime with an interval for new_regime (if they differ)
		if (index == m_dimension - 1 || is_changepoint_index_separator_index(index + 1)) {
			if (new_regimes[j] != removed_regimes[j]) {
				m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
				m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);
				m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
			}
		}
		else {
			unsigned int subs_regime_index = m_tau[index + 1].get_full_index_row(j);
			if (removed_regimes[j] != new_regimes[j]) {
				m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
				m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j], subs_regime_index);
				m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
			}
		}
        // calculate the trace index
        is_changepoint_index_separator_index(index); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index + 1) has just been run
        
        // if we are removing an unobserved regime and we don't add it back in, then the top regime must be removed
        if (removing_unobserved_regimes[j] && new_regimes[j] != m_regimes[j].size() - 1) {
            if (!m_unobserved_regimes[j].back()) {
                cerr << "regime to delete is not unobserved" << endl;
            }
            m_regimes[j].erase(m_regimes[j].end() - 1);
            unsigned int regime_to_remove = static_cast< unsigned int >(m_regimes[j].size());
            // for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
            for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
                m_regimes[j][regime].remove_unobserved_regime(regime_to_remove);
            }
            m_unobserved_regimes[j].erase(m_unobserved_regimes[j].end() - 1);
            m_number_of_unobserved_regimes[j]--;
        }
	}
}

void particle::resample_full_changepoint(const unsigned int & index, const vector< unsigned int > & new_regimes, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & regime_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes) {
    // calculate the regimes that we are going to lose, and set m_tau[index] to be associated with the new regimes.
    vector< unsigned int > removed_regimes(0);
    removed_regimes = m_tau[index].get_full_index();
    m_tau[index].set_full_index(new_regimes);
    size_t number_of_processes = new_regimes.size();
    
    for (unsigned int j = 0; j < number_of_processes; j++) {
        // if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
        bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
        if (adding_new_regime) {
            for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
                // increase the size of the transitions histogram for each of the regimes in this process.
                m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
            }
            m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(right_sufficient_statistics_reverse[j].size(), 0), m_number_of_traces, m_trace_index, 0));
            m_unobserved_regimes[j].push_back(true);
            m_number_of_unobserved_regimes[j]++;
        }
        
        unsigned int prev_regime_index;
        // calculate the regime index of the previous changepoint
        if (index == 0) {
            prev_regime_index = 0;
        }
        else {
            prev_regime_index = m_tau[index - 1].get_full_index_row(j);
        }
        
        // replace the interval for removed_regime with an interval for new_regime (if they differ)
        if (index == m_dimension - 1 || is_changepoint_index_separator_index(index + 1)) {
            if (index == m_dimension - 1) {
                // there is a separator changepoint at the beginning of the last trace, so we must be in the last trace.
                m_h_trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
            }
            else {
                // calculate the trace index
                is_changepoint_index_separator_index(index); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index + 1) has just been run
                m_h_trace_index = m_trace_index;
            }
            if (new_regimes[j] != removed_regimes[j]) {
                if (!is_changepoint_index_separator_index(index)) {
                    m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
                }
                m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);
                m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
            }
        }
        else {
            m_h_trace_index = m_trace_index;
            unsigned int subs_regime_index = m_tau[index + 1].get_full_index_row(j);
            if (removed_regimes[j] != new_regimes[j]) {
                m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
                m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j], subs_regime_index);
                if (!is_changepoint_index_separator_index(index)) {
                    m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
                }
            }
        }
        
        if (removed_regimes[j] != new_regimes[j]) {
            m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
            m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
            m_regimes[j][removed_regimes[j]].remove_observations(m_h_trace_index, right_number_of_observations_reverse[j]);
            m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics_reverse[j]);
            m_regimes[j][new_regimes[j]].set_log_likelihood(regime_log_likelihoods_with_right_sufficient_statistics_reverse[j][new_regimes[j]]);
            m_regimes[j][new_regimes[j]].add_observations(m_h_trace_index, right_number_of_observations_reverse[j]);
        }
        
        // if we are removing an unobserved regime and we don't add it back in, then the top regime must be removed
        if (removing_unobserved_regimes[j] && new_regimes[j] != m_regimes[j].size() - 1) {
            if (!m_unobserved_regimes[j].back()) {
                cerr << "regime to delete is not unobserved" << endl;
            }
            m_regimes[j].erase(m_regimes[j].end() - 1);
            unsigned int regime_to_remove = static_cast< unsigned int >(m_regimes[j].size());
            // for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
            for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
                m_regimes[j][regime].remove_unobserved_regime(regime_to_remove);
            }
            m_unobserved_regimes[j].erase(m_unobserved_regimes[j].end() - 1);
            m_number_of_unobserved_regimes[j]--;
        }
    }
}

void particle::alter_unobserved_regimes(const vector< int > & altering_unobserved_regimes, const size_t & number_of_processes) {
    for (unsigned int j = 0; j < number_of_processes; j++) {
        if (altering_unobserved_regimes[j] == 1) {
            // choose the new unobserved regime (the only restriction is that it can't be regime 0)
            unsigned int new_unobserved_regime = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_regimes[j].size()) + 1);
            size_t size_of_sufficient_statistics = m_regimes[j][0].get_size_of_sufficient_statistics();
            // insert the new unobserved regime into m_regimes[j]
            m_regimes[j].insert(m_regimes[j].begin() + new_unobserved_regime, regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size(), 0), vector< double >(size_of_sufficient_statistics, 0), m_number_of_traces, 0, 0));
            // set the likelihood to be 0 for this new regime
            m_regimes[j][new_unobserved_regime].set_log_likelihood(0);
            // for all other regimes, include an empty transition in m_right_transitions_histogram and, for each element in m_right_transitions, if the element is greater than or equal to the new_unobserved_regime then add 1 to it
            for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
                m_regimes[j][regime].insert_new_unobserved_regime(new_unobserved_regime);
            }
            for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
                m_tau[cp_index].increase_full_index_row_for_inserting_unobserved_regime_if_necessary(j, new_unobserved_regime);
            }
            m_unobserved_regimes[j].insert(m_unobserved_regimes[j].begin() + new_unobserved_regime, true);
            m_number_of_unobserved_regimes[j]++;
        }
        else if (altering_unobserved_regimes[j] == -1) {
            // choose which unobserved regime to remove using sub-linear method - keep guessing at unobserved regimes
            unsigned int regime_to_remove = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_regimes[j].size() - 1) + 1);
            bool unobserved = m_unobserved_regimes[j][regime_to_remove];
            while (!unobserved) {
                regime_to_remove = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_regimes[j].size() - 1) + 1);
                unobserved = m_unobserved_regimes[j][regime_to_remove];
            }
            m_regimes[j].erase(m_regimes[j].begin() + regime_to_remove);
            // for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
            for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
                m_regimes[j][regime].remove_unobserved_regime(regime_to_remove);
            }
            for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
                m_tau[cp_index].decrease_full_index_row_for_removing_unobserved_regime_if_necessary(j, regime_to_remove);
            }
            m_unobserved_regimes[j].erase(m_unobserved_regimes[j].begin() + regime_to_remove);
            m_number_of_unobserved_regimes[j]--;
        }
    }
}

//each changepoint now becomes the beginning of a binary
void particle::set_binary_marked_vectors(const size_t & number_of_processes, const vector< vector< double > > & log_likelihoods){
    for (size_t j = 0; j < number_of_processes; j++){
        m_binaries.push_back(vector< binary >(m_dimension + 2));
        m_binaries[j][0] = binary(log_likelihoods[j][0], -1);
        for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++){
            m_binaries[j][cp_index + 1] = binary(log_likelihoods[j][cp_index + 1], cp_index);
        }
        m_binaries[j][m_dimension + 1] = binary(-1e300, m_dimension);
    }
    m_intercept_changepoint.set_binary_index(vector< unsigned int >(number_of_processes, 0));
    for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++){
        m_tau[cp_index].set_binary_index(vector< unsigned int >(number_of_processes, cp_index + 1));
    }
}

void particle::set_all_binary_marked_vectors_equal_to_0_vectors(const size_t & number_of_processes) {
	m_intercept_changepoint.set_binary_index(vector< unsigned int >(number_of_processes, 0));
	for (int i = 0; i < m_dimension; i++) {
		m_tau[i].set_binary_index(vector< unsigned int >(number_of_processes, 0));
	}
}

void particle::add_end_binaries(const size_t & number_of_processes) {
	for (unsigned int foo = 0; foo < number_of_processes; foo++) {
		m_binaries[foo].push_back(binary(-1e300, m_dimension));
	}
}

// use the full changepoint indices in full_particle to set binary indices in this particle. Changepoints after the full_particle dimension will be set as new binaries
void particle::set_binary_indices_from_full_indices(particle & full_particle, const size_t & number_of_processes) {
    size_t full_dimension = full_particle.get_dimension();
    // we know that the full index and the binary index for the intercept changepoint must both be 0.
    m_binaries = vector< vector< binary > >(number_of_processes);
    for (unsigned int process = 0; process < number_of_processes; process++) {
        unsigned int binary_index = 0;
        // set the 0th binary to have likelihood = 0, left_index = -1
        m_binaries[process].push_back(binary(0, -1));
        if (full_dimension > 0) {
            if (full_particle.get_full_I_i_j(0, process) == 0) {
                // if the 0th changepoint has full index == 0, then we do not need to make a new binary
                m_tau[0].set_binary_index_row(process, binary_index);
            }
            else {
                // if the 0th changpoint is a change in regime then we need to increase the binary index and make a new binary.
                m_tau[0].set_binary_index_row(process, ++binary_index);
                // make a new binary with likelihood 0 and left index 0 (the likelihoods will be set later)
                m_binaries[process].push_back(binary(0, 0));
            }
        }
        
        for (unsigned int cp_index = 1; cp_index < full_dimension; cp_index++) {
            // Any regime change? If so, make a new binary
            if (full_particle.get_full_I_i_j(cp_index, process) == full_particle.get_full_I_i_j(cp_index - 1, process)) {
                m_tau[cp_index].set_binary_index_row(process, binary_index);
            }
            else {
                m_tau[cp_index].set_binary_index_row(process, ++binary_index);
                // set a new binary with log_likelihood 0 (it will be assigned later) and left index = cp_index
                m_binaries[process].push_back(binary(0, cp_index));
            }
        }
        
        if (m_dimension > full_dimension) {
            for (unsigned int cp_index = static_cast< unsigned int >(full_dimension); cp_index < m_dimension; cp_index++) {
                m_tau[cp_index].set_binary_index_row(process, ++binary_index);
                m_binaries[process].push_back(binary(m_tau[cp_index].get_log_likelihood(), cp_index));
            }
        }
        
        m_binaries[process].push_back(binary(-1e300, m_dimension));
    }
}

// set the log likelihoods for all the binaries
void particle::set_binary_log_likelihoods(const vector< vector< double > > & log_likelihoods, const size_t & number_of_processes) {
    for (unsigned int process = 0; process < number_of_processes; process++) {
        if (log_likelihoods[process].size() != m_binaries[process].size() - 1) { cout << "log_likelihoods and binaries sizes don't match" << endl; }
        for (unsigned int binary_index = 0; binary_index < log_likelihoods.size(); binary_index++) {
            m_binaries[process][binary_index].set_log_likelihood(log_likelihoods[process][binary_index]);
        }
    }
}

vector< unsigned long int > particle::calculate_and_get_changepoint_histogram(const unsigned long int & number_of_changepoint_bins){
	vector< unsigned long int > changepoints_vector(number_of_changepoint_bins, 0);
    unsigned long int old_bin_index = m_dimension; // set this so that the bin_index can't equal old_bin_index
	for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++){
		unsigned long int bin_index = (number_of_changepoint_bins * (m_tau[cp_index].get_position() - 1)) / m_end;
		if (bin_index != old_bin_index){ // do we already know that there's at least one cp in this bin?
			changepoints_vector[bin_index]++;
		}
        old_bin_index = bin_index;
	}
	return changepoints_vector;
}

vector< unsigned long int > particle::calculate_and_get_binary_changepoint_histogram(const unsigned long int & number_of_changepoint_bins, const size_t & number_of_processes){
    vector< unsigned long int > changepoints_vector(number_of_changepoint_bins * number_of_processes, 0);
    for (unsigned int process = 0; process < number_of_processes; process++) {
        size_t number_of_binaries = m_binaries[process].size() - 2;
        unsigned long int old_bin_index = number_of_binaries + 1; // set this so that the bin_index doesn't begin equalling the old_bin_index
        for (unsigned int binary_index = 1; binary_index < number_of_binaries + 1; binary_index++) {
            unsigned long int bin_index = (number_of_changepoint_bins * (m_tau[m_binaries[process][binary_index].get_left_index()].get_position() - 1)) / m_end + number_of_changepoint_bins * process;
            if (bin_index != old_bin_index) { // have we already recorded that there's at least one cp in this bin?
                changepoints_vector[bin_index]++;
            }
            old_bin_index = bin_index;
        }
    }
    return changepoints_vector;
}

vector< unsigned long int > particle::calculate_and_get_full_changepoint_histogram(const unsigned long int & number_of_changepoint_bins, const size_t & number_of_processes) {
	vector< unsigned long int > changepoints_vector(number_of_changepoint_bins * number_of_processes, 0);
	unsigned long int old_bin_index;
	unsigned long int bin_index;
	for (unsigned int process = 0; process < number_of_processes; process++) {
		// old_bin_index is a dummy to make sure that we don't keep adding changepoints from the same bin
		old_bin_index = -1;
        // is m_tau[0] an effective changepoint?
        if (0 < m_dimension && m_tau[0].get_full_index_row(process) != 0) {
            bin_index = (number_of_changepoint_bins * (m_tau[0].get_position() - 1)) / m_end + number_of_changepoint_bins * process;
            // have we already added a changepoint for this section?
            if (bin_index != old_bin_index) {
                changepoints_vector[bin_index]++;
            }
            old_bin_index = bin_index;
        }
		for (unsigned int cp_index = 1; cp_index < m_dimension; cp_index++) {
            if (m_tau[cp_index].get_full_index_row(process) != m_tau[cp_index - 1].get_full_index_row(process)) {
                bin_index = (number_of_changepoint_bins * (m_tau[cp_index].get_position() - 1)) / m_end + number_of_changepoint_bins * process;
                // have we already added a changepoint for this section?
                if (bin_index != old_bin_index) {
                    changepoints_vector[bin_index]++;
                }
                old_bin_index = bin_index;
            }
		}
	}
	return changepoints_vector;
}

void particle::calculate_and_add_similarity_matrices(vector< vector< double > > & similarity_matrix, const size_t & number_of_processes) {
    for (unsigned int proc = 0; proc < number_of_processes; proc++) {
        for (unsigned int reg_idx = 0; reg_idx < m_regimes[proc].size(); reg_idx++) {
            m_regimes[proc][reg_idx].add_to_similarity_matrix(similarity_matrix);
        }
    }
}

void particle::calculate_and_add_min_proportion_similarity_matrices(vector< vector< double > > & min_proportion_similarity_matrix, const size_t & number_of_processes, const vector< unsigned long int > & actual_number_of_observations) {
    // make two matrices, with 0 giving the sum over processes and regimes of I(n_{r,i,.}^{(s_0)} > 0, n_{r,i,.}^{(s_1)} > 0) n_{r,i,.}^{(s_0)} all divided by n_{s_0} and 1 giving the same but for the second trace. We will take the minimum after all computation.
    vector< vector< double > > min_prop_mat_0(min_proportion_similarity_matrix.size(), vector< double >(min_proportion_similarity_matrix[0].size(), 0));
    vector< vector< double > > min_prop_mat_1(min_proportion_similarity_matrix.size(), vector< double >(min_proportion_similarity_matrix[0].size(), 0));
    for (unsigned int proc = 0; proc < number_of_processes; proc++) {
        for (unsigned int reg_idx = 0; reg_idx < m_regimes[proc].size(); reg_idx++) {
            m_regimes[proc][reg_idx].add_to_min_proportion_similarity_matrices(min_prop_mat_0, min_prop_mat_1, actual_number_of_observations);
        }
    }
    // now work out which of these is the minimum for each pair of traces
    for (unsigned int ind_0 = 0; ind_0 < min_proportion_similarity_matrix.size(); ind_0++) {
        for (unsigned int ind_1 = 0; ind_1 < min_proportion_similarity_matrix.size(); ind_1++) {
            min_proportion_similarity_matrix[ind_0][ind_1] += min_prop_mat_0[ind_0][ind_1] < min_prop_mat_1[ind_0][ind_1] ? min_prop_mat_0[ind_0][ind_1] : min_prop_mat_1[ind_0][ind_1];
        }
    }
}

void particle::add_to_association_matrix(vector< vector< double > > & association_matrix, const unsigned int & process) {
    unsigned long int number_of_association_matrix_bins = association_matrix.size();
    // for each regime that affects this process, add implied associations to matrix
    for (unsigned int regime_index = 0; regime_index < m_regimes[process].size(); regime_index++) {
        vector< int > regime_right_cp_indices = m_regimes[process][regime_index].get_right_changepoint_indices();
        // work out the indices of the bins that are within this regime
        vector< unsigned long int > regime_bin_indices = vector< unsigned long int >(0);
        for (unsigned int regime_right_cp_indices_index = 0; regime_right_cp_indices_index < regime_right_cp_indices.size(); regime_right_cp_indices_index++) {
            // work out the left and right change point indices for this interval
            int right_cp_index = regime_right_cp_indices[regime_right_cp_indices_index];
            unsigned long int right_changepoint_position;
            if (right_cp_index == m_dimension) {
                right_changepoint_position = m_end;
            }
            else {
                right_changepoint_position = m_tau[right_cp_index].get_position();
            }
            unsigned long int left_changepoint_position;
            if (right_cp_index == 0) {
                left_changepoint_position = 0;
            }
            else {
                left_changepoint_position = m_tau[right_cp_index - 1].get_position();
            }
            // i is the first bin association bin point that lies within the range of this interval
            unsigned long int i, i_limit;
            if ((left_changepoint_position * number_of_association_matrix_bins) % m_end == 0) {
                i = (left_changepoint_position * number_of_association_matrix_bins) / m_end;
            }
            else {
                i = (left_changepoint_position * number_of_association_matrix_bins) / m_end + 1;
            }
            // now feed in valid values of i to regime_bin_indices
            if ((right_changepoint_position * number_of_association_matrix_bins) % m_end == 0) {
                i_limit = (right_changepoint_position * number_of_association_matrix_bins) / m_end;
            }
            else {
                i_limit = (right_changepoint_position * number_of_association_matrix_bins) / m_end + 1;
            }
            while (i < i_limit) {
                regime_bin_indices.push_back(i);
                i++;
            }
        }
        size_t number_of_bin_indices = regime_bin_indices.size();
        for (unsigned int bin_index = 0; bin_index < number_of_bin_indices; bin_index++) {
            for (unsigned int bin_index_1 = 0; bin_index_1 < number_of_bin_indices; bin_index_1++) {
                unsigned long int i_0 = regime_bin_indices[bin_index], i_1 = regime_bin_indices[bin_index_1];
                association_matrix[i_0][i_1] += 1;
            }
        }
    }
}

double particle::get_basic_log_posterior(){
	return m_log_likelihood + m_log_k_prior;
}

double particle::get_binary_log_posterior(){
    return m_log_likelihood + m_log_k_prior + m_log_binary_I_prior;
}

double particle::get_full_log_posterior(){
    return m_log_likelihood + m_log_k_prior + m_log_full_I_prior + m_log_regimes_prior + m_log_full_separators_prior;
}

double particle::get_full_log_SMC_posterior() {
    return m_log_likelihood + m_log_k_prior + m_log_full_I_prior + m_log_regimes_prior + m_log_full_separators_prior + m_log_arbitrary_extension_density;
}

// this function is called by a function manually calculating the log_likelihood for the process
vector< unsigned long int > particle::get_vector_of_binary_left_positions(const unsigned int process) {
    vector< unsigned long int > binary_left_positions;
    size_t number_of_left_positions = m_binaries[process].size();
    for (unsigned int i = 0; i < number_of_left_positions; i++) {
        int left_index = m_binaries[process][i].get_left_index();
        if (left_index == -1) {
            binary_left_positions.push_back(m_intercept_changepoint.get_position());
        }
        else if (left_index == m_dimension) {
            binary_left_positions.push_back(m_end + 1);
        }
        else {
            binary_left_positions.push_back(m_tau[left_index].get_position());
        }
    }
    return binary_left_positions;
}

// work out the left indices for all the changepoints that correspond to a left index for a binary for the given process
vector< int > particle::get_vector_of_binary_left_indices(const unsigned int process) {
    vector< int > binary_left_indices;
    size_t number_of_left_indices = m_binaries[process].size();
    for (unsigned int i = 0; i < number_of_left_indices; i++) {
        int left_index = m_binaries[process][i].get_left_index();
        binary_left_indices.push_back(left_index);
    }
    return binary_left_indices;
}

// get the left index of the binary containing m_tau[index - 1]. Often used to test if I_{h,j} == 1, as get_binary_left_index(j, h + 1) == h if this is true
int particle::get_binary_left_index(const size_t & process, const unsigned int & index){
    if (index == 0){
        return -1;
    } else {
        return m_binaries[process][m_tau[index - 1].get_binary_index_row(process)].get_left_index();
    }
}

// get the left index of the binary to the right of the one containing m_tau[index - 1]
int particle::get_binary_right_index(const size_t & process, const unsigned int & index){
    if (index == 0){
        return m_binaries[process][1].get_left_index();
    } else {
        return m_binaries[process][m_tau[index - 1].get_binary_index_row(process) + 1].get_left_index();
    }
}

// returns I_{i,j}
unsigned int particle::get_binary_I_i_j(const unsigned int & i, const unsigned int & j) { //i > -1
    if (i == m_binaries[j][m_tau[i].get_binary_index_row(j)].get_left_index()) { //if i is the beginning of a binary
        return 1;
    }
    else {
        return 0;
    }
}

// returns the index of the binary that contains m_tau[i] for process j
unsigned int particle::get_full_I_i_j(const unsigned int & i, const unsigned int & j) { //i > -1
    return m_tau[i].get_full_index_row(j);
}

// returns the log_likelihood of binary that contains m_tau[index - 1]
double particle::get_binary_log_likelihood(const unsigned int & j, const unsigned int & index){
    if (index == 0){
        return m_binaries[j][0].get_log_likelihood();
    } else {
        return m_binaries[j][m_tau[index - 1].get_binary_index_row(j)].get_log_likelihood();
    }
}

// sets the regimes based on the existing binaries and the sufficient statistics. sufficient statistics gives the vector vectors of cumulative data for each changepoint. e.g. sufficient_statistics[0][0] gives the cumulative sufficient statistics up to changepoint -1 for process 0.
void particle::set_full_marked_vectors(const size_t & number_of_processes, const vector< vector< vector< double > > > & sufficient_statistics, const vector< vector< double > > & number_of_observations) {
    m_log_regimes_prior = static_cast< double >(number_of_processes) * log(m_rho);
    m_log_full_separators_prior = 0;
    m_regimes = vector< vector< regime > >(0);
    m_unobserved_regimes = vector< vector< bool > >(0);
    // set the full index for each changepoint to be equal to the binary index for each changepoint
    m_intercept_changepoint.set_full_index_equal_to_binary_index();
    for (unsigned int i = 0; i < m_dimension; i++) {
        m_tau[i].set_full_index_equal_to_binary_index();
    }
    for (unsigned int j = 0; j < number_of_processes; j++) {
        // each binary corresponds to a regime
        m_regimes.push_back(vector< regime >(0));
        unsigned int trace_index = 0;
        for (unsigned int binary_index = 0; binary_index < m_binaries[j].size() - 2; binary_index++) {
            // calculate sufficient statistics
            vector< double > stats1 = sufficient_statistics[m_binaries[j][binary_index].get_left_index() + 1][j]; // +1 because sufficient stats contains the suffcient stats for cp_index -1
            vector< double > stats2 = sufficient_statistics[m_binaries[j][binary_index + 1].get_left_index() + 1][j];
            for (unsigned int proc = 0; proc < number_of_processes; proc++) {
                stats2[proc] -= stats1[proc];
            }
            // calculate number of observations
            double number_of_obs = number_of_observations[m_binaries[j][binary_index + 1].get_left_index() + 1][j] - number_of_observations[m_binaries[j][binary_index].get_left_index() + 1][j];
            // calculate the right changepoint indices for each regime
            vector< int > right_indices = vector< int >(0);
            for (int i = m_binaries[j][binary_index].get_left_index(); i < m_binaries[j][binary_index + 1].get_left_index(); i++) {
                right_indices.push_back(i + 1);
            }
            
            // calculate the right transitions and the right transitions histogram
            vector< unsigned int > transitions = vector< unsigned int >(0);
            vector< unsigned int > transitions_histogram = vector< unsigned int >(m_binaries[j].size() - 1, 0); // there is no slot for a non-regime at the end (can't transition to the non-regime at the end)
            for (int i = m_binaries[j][binary_index].get_left_index(); i < m_binaries[j][binary_index + 1].get_left_index() - 1; i++) {
                transitions.push_back(binary_index);
                transitions_histogram[binary_index]++;
            }
            if (is_changepoint_index_separator_index(m_binaries[j][binary_index + 1].get_left_index())) {
                transitions.push_back(-1);
            }
            else {
                transitions.push_back(binary_index + 1);
                transitions_histogram[binary_index + 1]++;
            }
            if (is_changepoint_index_separator_index(m_binaries[j][binary_index].get_left_index())) {
                trace_index++;
            }

            m_regimes[j].push_back(regime(right_indices, transitions, transitions_histogram, stats2, m_number_of_traces, trace_index, number_of_obs));
            // set the likelihood for the regime
            m_regimes[j][binary_index].set_log_likelihood(m_binaries[j][binary_index].get_log_likelihood());
        }
        // calculate sufficient statistics
        unsigned int binary_index = static_cast< unsigned int >(m_binaries[j].size() - 2);
        vector< double > stats1 = sufficient_statistics[m_binaries[j][binary_index].get_left_index() + 1][j]; // +1 because sufficient_stats contains the suffcient stats for cp_index -1
        vector< double > stats2 = sufficient_statistics[m_binaries[j][binary_index + 1].get_left_index() + 1][j];
        for (unsigned int i = 0; i < stats2.size(); i++) {
            stats2[i] -= stats1[i];
        }
        // calculate number of observations
        double number_of_obs = number_of_observations[m_binaries[j][binary_index + 1].get_left_index() + 1][j] - number_of_observations[m_binaries[j][binary_index].get_left_index() + 1][j];
        // calculate the left and right changepoint indices for each regime
        vector< int > right_indices = vector< int >(0);
        for (int i = m_binaries[j][binary_index].get_left_index(); i < m_binaries[j][binary_index + 1].get_left_index(); i++) {
            right_indices.push_back(i + 1);
        }
        // calculate the right transitions and the right transitions histogram
        vector< unsigned int > transitions = vector< unsigned int >(0);
        for (int i = m_binaries[j][binary_index].get_left_index(); i < m_binaries[j][binary_index + 1].get_left_index() - 1; i++) {
            transitions.push_back(binary_index);
        }
        // add a transition to nothing at the end.
        transitions.push_back(-1);
        vector< unsigned int > transitions_histogram = vector< unsigned int >(m_binaries[j].size() - 1, 0); // there is no slot for a non-regime at the end (can't transition to the non-regime at the end)
        transitions_histogram[binary_index] = m_binaries[j][binary_index + 1].get_left_index() - m_binaries[j][binary_index].get_left_index() - 1;
        // create the regime
        if (m_include_separator) {
            trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
        }
        m_regimes[j].push_back(regime(right_indices, transitions, transitions_histogram, stats2, m_number_of_traces, trace_index, number_of_obs));
        // set the log likelihood for the new regime
        m_regimes[j][binary_index].set_log_likelihood(m_binaries[j][binary_index].get_log_likelihood());
        
        m_log_regimes_prior += static_cast< double >(m_regimes[j].size() - 1) * log(1 - m_rho);
        m_log_full_separators_prior -= static_cast< double >(m_separator_indices.size()) * log(static_cast< double >(m_regimes[j].size()));
        // set the unobserved regimes indicator to be false for every regime
        m_unobserved_regimes.push_back(vector< bool >(m_regimes[j].size(), false));
    }
    // set the number of unobserved regimes for each process to be 0
    m_number_of_unobserved_regimes = vector< unsigned int >(number_of_processes, 0);
}

void particle::set_full_marked_vectors_without_binaries(const size_t & number_of_processes, const vector< vector< vector< double > > > & sufficient_statistics, const vector< vector< double > > & number_of_observations, const vector< vector< double > > & log_likelihoods) {
    m_log_regimes_prior = static_cast< double >(number_of_processes) * log(m_rho);
    m_log_full_separators_prior = 0;
    m_regimes = vector< vector< regime > >(0);
    m_unobserved_regimes = vector< vector< bool > >(0);
    // set the full index for each changepoint to be equal to the binary index for each changepoint
    m_intercept_changepoint.set_full_index(vector< unsigned int >(number_of_processes, 0));
    // add separator changepoints?
    for (unsigned int i = 0; i < m_dimension; i++) {
        m_tau[i].set_full_index(vector< unsigned int >(number_of_processes, i + 1));
    }
    for (unsigned int j = 0; j < number_of_processes; j++) {
        // each binary corresponds to a regime
        m_regimes.push_back(vector< regime >(0));
        unsigned int trace_index = 0;
        for (unsigned int regime_index = 0; regime_index < m_dimension; regime_index++) {
            // calculate sufficient statistics
            vector< double > stats1 = sufficient_statistics[regime_index][j];
            vector< double > stats2 = sufficient_statistics[regime_index + 1][j];
            for (unsigned int proc = 0; proc < stats1.size(); proc++) {
                stats2[proc] -= stats1[proc];
            }
            // calculate number of observations
            double number_of_obs = number_of_observations[regime_index + 1][j] - number_of_observations[regime_index][j];
            // calculate the right changepoint indices for each regime
            vector< int > right_indices = vector< int >(0);
            right_indices.push_back(regime_index);
            // calculate the right transitions and the right transitions histogram
            vector< unsigned int > transitions = vector< unsigned int >(0);
            vector< unsigned int > transitions_histogram = vector< unsigned int >(m_dimension + 1, 0); // there is no slot for a non-regime at the end (can't transition to the non-regime at the end)
            if (is_changepoint_index_separator_index(regime_index)) {
                transitions.push_back(-1);
            }
            else {
                transitions.push_back(regime_index + 1);
                transitions_histogram[regime_index + 1]++;
            }
            if (is_changepoint_index_separator_index(regime_index - 1)) {
                trace_index++;
            }
            
            m_regimes[j].push_back(regime(right_indices, transitions, transitions_histogram, stats2, m_number_of_traces, trace_index, number_of_obs));
            // set the likelihood for the regime
            m_regimes[j][regime_index].set_log_likelihood(log_likelihoods[regime_index + 1][j]);
        }
        // calculate sufficient statistics
        unsigned int regime_index = m_dimension;
        vector< double > stats1 = sufficient_statistics[regime_index][j];
        vector< double > stats2 = sufficient_statistics[regime_index + 1][j];
        for (unsigned int i = 0; i < stats2.size(); i++) {
            stats2[i] -= stats1[i];
        }
        // calculate number of observations
        double number_of_obs = number_of_observations[regime_index + 1][j] - number_of_observations[regime_index][j];
        // calculate the left and right changepoint indices for each regime
        vector< int > right_indices = vector< int >(0);
        right_indices.push_back(regime_index);
        // calculate the right transitions and the right transitions histogram
        vector< unsigned int > transitions = vector< unsigned int >(0);
        vector< unsigned int > transitions_histogram = vector< unsigned int >(m_dimension + 1, 0);
        // add a transition to nothing at the end.
        transitions.push_back(-1);
        // create the regime
        if (m_include_separator) {
            trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
        }
        m_regimes[j].push_back(regime(right_indices, transitions, transitions_histogram, stats2, m_number_of_traces, trace_index, number_of_obs));
        // set the log likelihood for the new regime
        m_regimes[j][regime_index].set_log_likelihood(log_likelihoods[regime_index + 1][j]); // IS THIS RIGHT?
        
        m_log_regimes_prior += static_cast< double >(m_regimes[j].size() - 1) * log(1 - m_rho);
        m_log_full_separators_prior -= static_cast< double >(m_separator_indices.size()) * log(static_cast< double >(m_regimes[j].size()));
        // set the unobserved regimes indicator to be false for every regime
        m_unobserved_regimes.push_back(vector< bool >(m_regimes[j].size(), false));
    }
    // set the number of unobserved regimes for each process to be 0
    m_number_of_unobserved_regimes = vector< unsigned int >(number_of_processes, 0);
}

void particle::set_all_full_marked_vectors_equal_to_binary_marked_vectors(const size_t & number_of_processes) {
    m_intercept_changepoint.set_full_index(vector< unsigned int >(number_of_processes, 0));
    // add separator changepoints?
    for (unsigned int i = 0; i < m_dimension; i++) {
        m_tau[i].set_full_index(vector< unsigned int >(number_of_processes, i + 1));
    }
}

void particle::set_all_regimes_to_be_observed(const size_t & number_of_processes) {
    m_unobserved_regimes = vector< vector< bool > >(0);
    m_number_of_unobserved_regimes = vector< unsigned int >(0);
    for (unsigned int foo = 0; foo < number_of_processes; foo++) {
        m_unobserved_regimes.push_back(vector< bool >(m_regimes[foo].size(), false));
        m_number_of_unobserved_regimes.push_back(0);
    }
}

void particle::add_new_regime(const unsigned int & process, vector< int > & right_changepoint_indices, vector< unsigned int > & right_transitions, vector< unsigned int > & right_transitions_histogram, vector< double > & sufficient_statistics, double & number_of_observations, double & log_likelihood, size_t & number_of_traces, const bool & new_trace, const unsigned int & trace_index, const unsigned int & previous_regime) {
    unsigned int new_regime_index = static_cast< unsigned int >(m_regimes[process].size());
    // increase the number of slots in the transitions_histogram for each other regime, adding a new transition from the previous regime unless this is the beginning of a new trace
    for (unsigned int bar = 0; bar < new_regime_index; bar++) {
        m_regimes[process][bar].add_regime_to_transitions_histogram(new_regime_index);
        if (bar == previous_regime && !new_trace) {
            m_regimes[process][bar].alter_last_transition(new_regime_index);
        }
    }
    m_regimes[process].push_back(regime(right_changepoint_indices, right_transitions, right_transitions_histogram, sufficient_statistics, number_of_traces, trace_index, number_of_observations));
    m_regimes[process].back().set_log_likelihood(log_likelihood);
    // alter the full_indices for the changepoints
    for (int foo = right_changepoint_indices[0]; foo <= right_changepoint_indices.back(); foo++) {
        if (0 < foo) {
            m_tau[foo - 1].set_full_index_row(process, new_regime_index);
        }
    }
    
    m_regimes[process].back().set_transitions_out(m_regimes[process].back().get_number_of_right_transitions());
    // the new regime will be set to be observed later
}

void particle::add_binary_to_regime(const unsigned int & process, const unsigned int & regime_index, const vector< int > & right_changepoint_indices, const vector< unsigned int > & right_transitions, const vector< unsigned int > & right_transitions_histogram, const vector< double > & extra_sufficient_statistics, const unsigned int & extra_number_of_observations, const double & log_likelihood_with_right_sufficient_statistics, const size_t & number_of_traces, const bool & new_trace, const unsigned int & trace_index, const unsigned int & previous_regime) {
    m_regimes[process][regime_index].append_right_changepoint_indices(right_changepoint_indices);
    m_regimes[process][regime_index].append_right_transitions(right_transitions);
    m_regimes[process][regime_index].add_right_transitions_histogram(right_transitions_histogram);
    m_regimes[process][regime_index].add_sufficient_statistics(extra_sufficient_statistics);
    m_regimes[process][regime_index].add_observations(trace_index, extra_number_of_observations);
    m_regimes[process][regime_index].set_log_likelihood(log_likelihood_with_right_sufficient_statistics);
    
    // alter the full_indices of the changepoints
    for (int foo = right_changepoint_indices[0]; foo <= right_changepoint_indices.back(); foo++) {
        if (0 < foo) {
            m_tau[foo - 1].set_full_index_row(process, regime_index);
        }
    }
    
    // if this is not a new trace, need to add a transition to the end of the previous regime
    if (!new_trace && previous_regime != regime_index) {
        m_regimes[process][previous_regime].alter_last_transition(regime_index);
    }
    else if (!new_trace && previous_regime == regime_index) {
        m_regimes[process][previous_regime].alter_right_transition(right_changepoint_indices[0] - 1, regime_index);
    }
}

// calculate and get the number of observed regimes
unsigned int particle::get_number_of_observed_regimes(const unsigned int & process) {
    unsigned int num_obs_regimes = 0;
    for (const auto obs_bool:m_unobserved_regimes[process]) {
        if (!obs_bool) {
            num_obs_regimes++;
        }
    }
    // check this matches number of regimes - number of unobserved regimes
    unsigned int x = static_cast< unsigned int >(m_regimes[process].size()) - m_number_of_unobserved_regimes[process];
    if (x != num_obs_regimes) {
        cerr << "number of unobserved regimes don't match!" << endl;
    }
    return num_obs_regimes;
}

bool particle::removing_full_changepoint_leaves_highest_regime_unobserved(const unsigned int & process, const unsigned int & regime) {
    return (m_regimes[process][regime].get_number_of_left_indices() == 1) && (regime > 0) && (regime == m_regimes[process].size() - 1);
}

bool particle::removing_full_changepoint_leaves_regime_unobserved(const unsigned int & process, const unsigned int & regime) {
	return (m_regimes[process][regime].get_number_of_left_indices() == 1) && (regime > 0);
}

// are there any unobserved regimes for this process?
vector< bool > particle::any_unobserved_regimes() {
    size_t number_of_processes = m_number_of_unobserved_regimes.size();
    vector< bool > any_unobserved_regimes(number_of_processes, 0);
    for (unsigned int j = 0; j < number_of_processes; j++) {
        any_unobserved_regimes[j] = m_number_of_unobserved_regimes[j] > 0;
    }
    return any_unobserved_regimes;
}

// find the regime index for the regime which affects the changepoint prior to cp_index
unsigned int particle::get_previous_regime(const int & cp_index, const unsigned int & process) {
    if (cp_index == 0) {
        return 0;
    }
    else {
        return m_tau[cp_index - 1].get_full_index_row(process);
    }
}

// obtain the vector of sufficient statistics for the regime that affects process with index regime_index
vector< double > particle::get_sufficient_statistics(const unsigned int & process, const unsigned int & regime_index) {
	return m_regimes[process][regime_index].get_sufficient_statistics();
}

// recover the log likelihood for regime regime_index for process
double particle::get_regime_log_likelihood(const unsigned int & process, const unsigned int & regime_index) {
	return m_regimes[process][regime_index].get_log_likelihood();
}

bool particle::deleting_left_changepoint_removes_regime(const unsigned int & process, const unsigned int & regime_index) {
    return m_regimes[process][regime_index].get_number_of_left_indices() > 1;
}

void particle::increase_separator_indices_greater_or_equal_to_index(const unsigned int & index) {
    unsigned int i = static_cast< unsigned int >(m_separator_indices.size()) - 1;
    bool greater_or_equal = m_separator_indices[i] >= index;
    bool can_go_lower = true;
    while (greater_or_equal && can_go_lower) {
        m_separator_indices[i]++;
        i--;
        can_go_lower = i != numeric_limits< unsigned int >::max();
        if (can_go_lower) {
            greater_or_equal = m_separator_indices[i] >= index;
        }
    }
}

void particle::decrease_separator_indices_greater_than_index(const unsigned int & index) {
    unsigned int i = static_cast< unsigned int >(m_separator_indices.size()) - 1;
    bool greater_or_equal = m_separator_indices[i] > index;
    bool can_go_lower = true;
    while (greater_or_equal && can_go_lower) {
        m_separator_indices[i]--;
        i--;
        can_go_lower = i != numeric_limits< unsigned int >::max();
        if (can_go_lower) {
            greater_or_equal = m_separator_indices[i] > index;
        }
    }
}

void particle::check_unobserved_regimes(const unsigned int process) {
    unsigned int no_unob_regs = 0;
    for (unsigned int i = 0; i < m_unobserved_regimes[process].size(); i++) {
        if (m_unobserved_regimes[process][i]) {
            no_unob_regs++;
        }
    }
    if (no_unob_regs != m_number_of_unobserved_regimes[process]) {
        cerr << "don't match" << endl;
    }
}

void particle::check_separator_changepoints() {
    if (m_separator_indices.size() != m_separators.size()) {
        cerr << "sizes don't match" << endl;
    }
    for (unsigned int i = 0; i < m_separator_indices.size(); i++) {
        if (m_separators[i] != m_tau[m_separator_indices[i]].get_position()) {
            cerr << "m_separators are wrong" << endl;
        }
    }
}

void particle::check_observations_in_traces(const unsigned long int & time) {
    double total_number_of_observations = 0;
    for (unsigned int process = 0; process < m_regimes.size(); process++) {
        for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
            for (unsigned int trace = 0; trace < m_number_of_traces; trace++) {
                if (m_regimes[process][regime].get_number_of_observations(trace) < -0.0001) {
                    cerr << "negative number of observations" << endl;
                }
            }
            double number_of_observations = 0;
            for (unsigned int trace = 0; trace < m_number_of_traces; trace++) {
                number_of_observations += m_regimes[process][regime].get_number_of_observations(trace);
            }
            total_number_of_observations += number_of_observations;
            double number_of_sufficient_statistics = 0;
            for (unsigned int i = 0; i < m_regimes[process][regime].get_sufficient_statistics().size(); i++) {
                number_of_sufficient_statistics += m_regimes[process][regime].get_sufficient_statistics()[i];
            }
            if (abs(number_of_sufficient_statistics - number_of_observations) > 0.0001) {
                cerr << "number of observations doesn't match" << endl;
            }
        }
    }
    /*if (abs(time - 1 - total_number_of_observations - static_cast< double >(m_number_of_traces - 1)) > 0.0001) {
        cerr << "wrong number of observations" << endl;
        cerr << "total number of observations: " << total_number_of_observations << ", time: " << time << ", m_trace_index: " << m_trace_index << endl;
    }*/
}

void particle::check_regimes_and_right_changepoints() {
    for (unsigned int proc = 0; proc < m_regimes.size(); proc++) {
        // check if regime 0 has changepoint 0 (possibly at end) as a right changepoint index
        if (m_regimes[proc][0].get_right_changepoint_indices()[0] != 0) {
            cerr << "regime 0 doesn't have changepoint 0 as a right changepoint index" << endl;
        }
        // go though each of the right changepoint indices in each regime and see if it matches the full indices of the changepoints
        for (unsigned int reg = 0; reg < m_regimes[proc].size(); reg++) {
            vector< int > right_changepoint_indices = m_regimes[proc][reg].get_right_changepoint_indices();
            vector< unsigned int > right_transitions = m_regimes[proc][reg].get_right_transitions();
            for (unsigned int i = 0; i < right_changepoint_indices.size(); i++) {
                if (right_changepoint_indices[i] != 0) {
                    if (reg != m_tau[right_changepoint_indices[i] - 1].get_full_index_row(proc)) {
                        cerr << "regimes don't match" << endl;
                    }
                }
                if (right_changepoint_indices[i] != m_dimension) {
                    if (right_transitions[i] != numeric_limits< unsigned int >::max()) {
                        if (right_transitions[i] != m_tau[right_changepoint_indices[i]].get_full_index_row(proc)) {
                            cerr << "regimes don't match" << endl;
                        }
                    }
                }
            }
            // also check that right transitions and right transitions histogram match up
            vector< unsigned int > right_transitions_histogram = m_regimes[proc][reg].get_right_transitions_histogram();
            for (unsigned int i = 0; i < right_transitions.size(); i++) {
                if (right_transitions[i] < numeric_limits< unsigned int >::max()) {
                    right_transitions_histogram[right_transitions[i]]--;
                }
            }
            for (unsigned int i = 0; i < right_transitions_histogram.size(); i++) {
                if (right_transitions_histogram[i] != 0) {
                    cerr << "transitions don't match" << endl;
                }
            }
        }
    }
}

void particle::check_full_log_posterior(const bool & always_new_regime = false) {
    unsigned int number_of_processes = static_cast< unsigned int >(m_regimes.size());
    double log_k_prior = 0;
    log_k_prior = calculate_and_get_log_k_prior();
    if (0.01 < abs(log_k_prior - m_log_k_prior)) {
        cout << "k prior difference " << log_k_prior - m_log_k_prior << endl;
    }
    double full_log_I_prior;
    if (always_new_regime) {
        full_log_I_prior = 0;
    }
    else {
        full_log_I_prior = calculate_and_get_log_full_I_prior(number_of_processes);
    }
    if (0.01 < abs(full_log_I_prior - m_log_full_I_prior)) {
        cout << "full log I prior difference " << full_log_I_prior - m_log_full_I_prior << endl;
    }
    
    double full_log_regime_prior;
    full_log_regime_prior = calculate_and_get_log_regimes_prior(number_of_processes);
    if (0.01 < abs(full_log_regime_prior - m_log_regimes_prior)) {
        cout << "full log regime prior difference " << full_log_regime_prior - m_log_regimes_prior << endl;
    }
    double full_log_separators_prior;
    full_log_separators_prior = calculate_and_get_log_full_separators_prior(number_of_processes);
    if (0.01 < abs(full_log_separators_prior - m_log_full_separators_prior)) {
        cout << "full log separators difference " << full_log_separators_prior - m_log_full_separators_prior << endl;
    }
    double full_log_likelihood;
    full_log_likelihood = calculate_and_get_full_log_likelihood(number_of_processes);
    if (0.01 < abs(full_log_likelihood - m_log_likelihood)) {
        cout << "full log likelihood difference " << full_log_likelihood - m_log_likelihood << endl;
    }
}

void particle::check_transitions_out() {
    for (unsigned int process = 0; process < m_regimes.size(); process++) {
        for (unsigned int reg = 0; reg < m_regimes[process].size(); reg++) {
            vector< unsigned int > right_transitions = m_regimes[process][reg].get_right_transitions_histogram();
            unsigned int trans_out = 0;
            for (unsigned int i = 0; i < right_transitions.size(); i++) {
                trans_out += right_transitions[i];
            }
            if (trans_out != m_regimes[process][reg].get_transitions_out()) {
                cout << "difference" << endl;
            }
        }
    }
}

//smc-only functions
double particle::calculate_and_get_add_cp_proposal_ratio(const unsigned long int & available_cp_positions, const unsigned long int & total_cp_positions) {
    double proposal_ratio = log(static_cast< double >(available_cp_positions)) - log(static_cast< double >(total_cp_positions - available_cp_positions + 1));
    if (total_cp_positions == available_cp_positions) { // i.e. dimension == 0
        proposal_ratio -= log(4); // prob 1 of proposing increase, prob 1/4 of proposing decrease, (1/4)/1 = 1/4
    }
    else if (available_cp_positions == 1) { // i.e. there is a single spot left for adding a changepoint
        proposal_ratio += log(4); // prob 1/4 of proposing increase, prob 1 of proposing decrease, 1/(1/4) = 4
    }
    return proposal_ratio;
}

double particle::calculate_and_get_remove_cp_proposal_ratio(const unsigned long int & available_cp_positions, const unsigned long int & total_cp_positions) {
    return log(static_cast< double >(total_cp_positions - available_cp_positions)) - log(available_cp_positions + 1);
}

void particle::extend_regime(const unsigned int & process, const unsigned int & previous_regime, const vector< double > & new_interval_sufficient_statistics, const double & previous_regime_new_log_likelihood, const unsigned int & trace_index, const double & number_of_observations) {
    m_regimes[process][previous_regime].add_sufficient_statistics(new_interval_sufficient_statistics);
	m_log_likelihood += previous_regime_new_log_likelihood - m_regimes[process][previous_regime].get_log_likelihood();
    m_regimes[process][previous_regime].set_log_likelihood(previous_regime_new_log_likelihood);
    m_regimes[process][previous_regime].add_observations(trace_index, number_of_observations);
}

void particle::add_full_separator_changepoint(changepoint & separator_changepoint, const vector< unsigned int > & new_regimes, const vector< vector< double > > & sufficient_statistics, const vector< vector< double > > & log_likelihoods, const vector< double > & number_of_observations) {
	m_log_k_prior += calculate_and_get_add_cp_k_prior_ratio();
	separator_changepoint.set_full_index(new_regimes);
	// insert new_changepoint into m_tau
    m_tau.push_back(separator_changepoint);
    m_dimension++; //increase dimension
    size_t number_of_processes = new_regimes.size();
    // make sure the m_include_separator is true
    m_include_separator = true;
	// add this separator to the list of separator indices
	m_separators.push_back(separator_changepoint.get_position());
	m_separator_indices.push_back(m_dimension - 1);
    m_number_of_traces++;
    // add an extra trace to each regime (for the number of observations in each regime)
    for (unsigned int proc = 0; proc < number_of_processes; proc++) {
        for (unsigned int reg = 0; reg < m_regimes[proc].size(); reg++) {
            m_regimes[proc][reg].add_new_trace();
        }
    }
	double s = static_cast<double>(m_separator_indices.size());
	for (unsigned int j = 0; j < number_of_processes; j++) {
        double number_of_regimes = static_cast< double >(m_regimes[j].size());
        m_log_full_separators_prior -= log(number_of_regimes);
		// if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
		bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
		if (adding_new_regime) {
			for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
				// increase the size of the transitions histogram for each of the regimes in this process.
				m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
			}
			for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
				m_log_full_I_prior += gsl_sf_lngamma((number_of_regimes + 1) * m_dirichlet_alpha) + gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + static_cast< double >(m_regimes[j][regime].get_number_of_right_transitions())) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma((number_of_regimes + 1) * m_dirichlet_alpha + static_cast< double >(m_regimes[j][regime].get_number_of_right_transitions()));
			}
			m_log_regimes_prior += log(1.0 - m_rho);
			m_log_full_separators_prior += s * log(number_of_regimes) - s * log(number_of_regimes + 1);
			// insert new regime into m_regimes[j]. We don't include the right changepoint position, right transition, sufficient statistics or number_of_observations here as it will be added later
			// these will be converted to the correct values later
			m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(sufficient_statistics[j].size(), 0), m_number_of_traces, m_trace_index, 0));
			// these unobserved regime values will be corrected later
			m_unobserved_regimes[j].push_back(true);
			m_number_of_unobserved_regimes[j]++;
		}
	}
    
    m_trace_index = static_cast< unsigned int >(m_separators.size());
    for (unsigned int j = 0; j < number_of_processes; j++) {
        // insert the interval
        m_regimes[j][new_regimes[j]].add_interval(m_dimension, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);
		m_log_likelihood += log_likelihoods[j][new_regimes[j]] - m_regimes[j][new_regimes[j]].get_log_likelihood();
        // assign right sufficient statistics and observation numbers to the new regime. Set the log likelihood for each regime
        m_regimes[j][new_regimes[j]].add_sufficient_statistics(sufficient_statistics[j]);
        m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods[j][new_regimes[j]]);
        m_regimes[j][new_regimes[j]].add_observations(m_trace_index, number_of_observations[j]);
    }
}

// add a changepoint that was guaranteed to be accepted to the particle. It will be added to the end of m_tau and new sufficient statistics, observations and likelihood will be added
void particle::add_guaranteed_full_changepoint(changepoint & adding_changepoint, const vector< unsigned int > & new_regimes, const vector< vector< double > > & right_sufficient_statistics, const vector< vector< double > > & log_likelihoods_with_right_sufficient_statistics, const vector< double > & right_number_of_observations) {
    m_log_k_prior += calculate_and_get_add_cp_k_prior_ratio();
	// insert new_changepoint into m_tau
	m_tau.push_back(adding_changepoint);
	m_dimension++; //increase dimension
	size_t number_of_processes = new_regimes.size();
	// assign the regime indices for the new changepoint.
	m_tau.back().set_full_index(new_regimes);
	for (unsigned int j = 0; j < number_of_processes; j++) {
		unsigned int prev_regime_index;
		// calculate the regime index of the previous changepoint
		if (m_dimension == 1) {
			prev_regime_index = 0;
		}
		else {
			prev_regime_index = m_tau[m_dimension - 2].get_full_index_row(j);
		}
		// if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
		bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
		if (adding_new_regime) {
			for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
				// increase the size of the transitions histogram for each of the regimes in this process.
				m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
			}
			// insert new regime into m_regimes[j]. We don't include the right changepoint position, right transition, sufficient statistics or number_of_observations here as it will be added later
			// these will be converted to the correct values later
			m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(right_sufficient_statistics[j].size(), 0), m_number_of_traces, m_trace_index, 0));
			// these unobserved regime values will be corrected later
			m_unobserved_regimes[j].push_back(true);
			m_number_of_unobserved_regimes[j]++;
		}
		// insert the interval and transition out of the new_regime
		// replace the right transition from previous regime to subsequent regime with a transition to new_regime (unless new_regime is equal to subsequent_regime), adding a transition out of previous regime if the new changepoint is inserted at the end of m_tau.
		m_regimes[j][prev_regime_index].alter_right_transition(m_dimension - 1, new_regimes[j]);
		m_regimes[j][new_regimes[j]].add_interval(m_dimension, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);
		m_trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
		
		// assign right sufficient statistics and observation numbers to the new regime. Calculate the log likelihood for each regime
        m_log_likelihood += log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]] - m_regimes[j][new_regimes[j]].get_log_likelihood();
		m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics[j]);
		m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]]);
		m_regimes[j][new_regimes[j]].add_observations(m_trace_index, right_number_of_observations[j]);
	}
}

void particle::increase_log_separator_prior(const double & full_log_prior, const bool & adding_new_regime, const unsigned int & process) {
    // if we are not adding a new regime then full_log_prior = -log(number_of_regimes)
    //m_log_full_separators_prior += full_log_prior;
    if (adding_new_regime) {
        // if we are adding a new regime then full_log_prior = log_full_I_prior_ratio + log(1 - rho) + s * log(number_of_regimes) - (s + 1) * log(number_of_regimes + 1)
        //m_log_full_separators_prior -= log(1.0 - m_rho);
        m_log_regimes_prior += log(1 - m_rho);
        double number_of_regimes = static_cast< double >(m_regimes[process].size());
        double s = static_cast<double>(m_separator_indices.size());
        double log_full_I_prior_ratio = full_log_prior - log(1.0 - m_rho) - s * log(number_of_regimes) + (s + 1) * log(number_of_regimes + 1);
        m_log_full_I_prior += log_full_I_prior_ratio;
        //m_log_full_separators_prior -= log_full_I_prior_ratio;
    }
}

// can either pass uniform RVs to choose the number of regimes for each process or pass in the number of unobserved regimes for each process
void particle::resample_number_of_unobserved_regimes(vector< double > uniform_rvs, vector< unsigned int > new_numbers_of_unobserved_regimes = vector< unsigned int >(0)) { 
	unsigned int number_of_processes = static_cast< unsigned int >(m_regimes.size());
	for (unsigned int proc = 0; proc < number_of_processes; proc++) {
		unsigned int new_number_of_unobserved_regimes = 0;
		if (0 < new_numbers_of_unobserved_regimes.size()) { // is new_numbers_of_unobserved_regimes not an empty vector?
			new_number_of_unobserved_regimes = new_numbers_of_unobserved_regimes[proc];
		}
		else {
			// use the uniform RV to choose the number of unobserved regimes for this process
			if (m_rho < uniform_rvs[proc]) {
				new_number_of_unobserved_regimes = 1;
				double temp = m_rho * (1 - m_rho);
				double temp_sum = m_rho + temp;
				while (temp_sum < uniform_rvs[proc]) {
					new_number_of_unobserved_regimes++;
					temp *= (1 - m_rho);
					temp_sum += temp;
				}
			}
		}
		if (m_number_of_unobserved_regimes[proc] < new_number_of_unobserved_regimes) { // need to add more unobserved regimes
			for (unsigned int new_reg = 0; new_reg < new_number_of_unobserved_regimes - m_number_of_unobserved_regimes[proc]; new_reg++) {
				// alter the priors: full_I, separators, regimes
                double number_of_regimes = static_cast< double >(m_regimes[proc].size());
				double log_full_I_prior_ratio = calculate_and_get_add_unobserved_regimes_full_I_prior_ratio(proc) - static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
                m_log_full_I_prior += log_full_I_prior_ratio;
                m_log_full_separators_prior += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
                m_log_regimes_prior += log(1 - m_rho);
				// add the new regime
				unsigned int new_regime = static_cast< unsigned int >(m_regimes[proc].size());
				for (unsigned int i = 0; i < m_regimes[proc].size(); i++) {
					// increase the size of the transitions histogram for each of the regimes in this process.
					m_regimes[proc][i].add_regime_to_transitions_histogram(new_regime);
				}
				// calculate the size of the sufficient statistics vector for process proc
				size_t suff_stats_size = m_regimes[proc][0].get_size_of_sufficient_statistics();
				// insert new regime into m_regimes[proc]
				m_regimes[proc].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[proc].size() + 1, 0), vector< double >(suff_stats_size, 0), m_number_of_traces, m_trace_index, 0));
				m_unobserved_regimes[proc].push_back(true);
				m_number_of_unobserved_regimes[proc]++;
			}
		}
		else if (new_number_of_unobserved_regimes < m_number_of_unobserved_regimes[proc]) { // removing unobserved regimes
			// find the indices of the unobserved regimes
			vector < unsigned int > unobserved_regime_indices = vector< unsigned int >(m_number_of_unobserved_regimes[proc]);
			unsigned int uo_idx = 0;
			for (unsigned int id = 0; id < m_unobserved_regimes[proc].size(); id++) {
				if (m_unobserved_regimes[proc][id]) {
					unobserved_regime_indices[uo_idx] = id;
					uo_idx++;
				}
			}
			for (int del_reg = new_number_of_unobserved_regimes - m_number_of_unobserved_regimes[proc]; del_reg < 0; del_reg++) {
                // alter the priors: full_I, separators, regimes
                double number_of_regimes = static_cast< double >(m_regimes[proc].size());
                double log_full_I_prior_ratio = calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio(proc) - static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes - 1));
                m_log_full_I_prior += log_full_I_prior_ratio;
                m_log_full_separators_prior += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes - 1));
                m_log_regimes_prior -= log(1 - m_rho);
				// delete the highest unobserved regime (so we don't have to keep changing the unobserved_regime_indices)
				m_regimes[proc].erase(m_regimes[proc].begin() + unobserved_regime_indices.back());
				unsigned int regime_to_remove = static_cast< unsigned int >(unobserved_regime_indices.back());
				unobserved_regime_indices.erase(unobserved_regime_indices.end() - 1);
				// for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
				for (unsigned int regime = 0; regime < m_regimes[proc].size(); regime++) {
					m_regimes[proc][regime].remove_unobserved_regime(regime_to_remove);
				}
                for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
                    m_tau[cp_index].decrease_full_index_row_for_removing_unobserved_regime_if_necessary(proc, regime_to_remove);
                }
				m_unobserved_regimes[proc].erase(m_unobserved_regimes[proc].begin() + regime_to_remove);
				m_number_of_unobserved_regimes[proc]--;
			}
		}
	}
}

void particle::permute_process_regimes(const unsigned int & process, const vector< unsigned int > & permutation_vector) {
	// permute the regimes stored in m_regimes[process] and the boolean values of m_unobserved_regimes
	vector< regime > regimes_copy = m_regimes[process];
	vector< bool > unobserved_regimes_copy = m_unobserved_regimes[process];
	for (unsigned int i = 0; i < m_regimes[process].size(); i++) {
		m_regimes[process][permutation_vector[i]] = regimes_copy[i];
		m_unobserved_regimes[process][permutation_vector[i]] = unobserved_regimes_copy[i];
	}
    for (unsigned int i = 0; i < m_regimes[process].size(); i++) {
        // permute the transitions histogram and right_transitions for this regime
        m_regimes[process][i].permute_transitions_histogram(permutation_vector);
        m_regimes[process][i].permute_right_transitions(permutation_vector);
    }
	// alter the full regime index for process for any changepoints (don't need to do the intercept changepoint because we know that the full index will be [0, ..., 0] for it.
	for (unsigned int j = 0; j < m_dimension; j++) {
		unsigned int current_reg = m_tau[j].get_full_index_row(process);
		m_tau[j].set_full_index_row(process, permutation_vector[current_reg]);
	}
}

bool particle::is_regime_observed_before(const unsigned int & cp_index, const unsigned int & process) { // is the regime for changepoint cp_index for the process observed before changepoint cp_index?
	unsigned int regime = m_tau[cp_index].get_full_index_row(process);
	if (regime == 0) {
		return true;
	}
	bool found = false;
	unsigned int j = 0;
	while (!found && j < cp_index) {
		found = m_tau[j].get_full_index_row(process) == regime;
		j++;
	}
	return found;
}

int particle::find_changepoint_with_same_regime(const unsigned int & cp_index, const unsigned int & process) { // given that the regime for the changepoint cp_index is observed before changepoint cp_index, when is it first observed?
	unsigned int regime = m_tau[cp_index].get_full_index_row(process);
	if (regime == 0) {
		return -1;
	}
	bool found = false;
	int j = 0;
	while (!found && j < cp_index) {
		found = m_tau[cp_index].get_full_index_row(process) == regime;
		j++;
	}
	j--;
	return j;
}

void particle::add_to_association_matrix_up_to_time(vector< vector< double > > & association_matrix, const unsigned int & process, const unsigned long int & time) {
    unsigned long int number_of_association_matrix_bins = association_matrix.size();
    // for each regime that affects this process, add implied associations to matrix
    for (unsigned int regime_index = 0; regime_index < m_regimes[process].size(); regime_index++) {
        vector< int > regime_right_cp_indices = m_regimes[process][regime_index].get_right_changepoint_indices();
        // work out the indices of the bins that are within this regime
        vector< unsigned long int > regime_bin_indices = vector< unsigned long int >(0);
        for (unsigned int regime_right_cp_indices_index = 0; regime_right_cp_indices_index < regime_right_cp_indices.size(); regime_right_cp_indices_index++) {
            // work out the left and right change point indices for this interval
            int right_cp_index = regime_right_cp_indices[regime_right_cp_indices_index];
            unsigned long int right_changepoint_position;
            if (right_cp_index == m_dimension) {
                right_changepoint_position = time;
            }
            else {
                right_changepoint_position = m_tau[right_cp_index].get_position();
            }
            unsigned long int left_changepoint_position;
            if (right_cp_index == 0) {
                left_changepoint_position = 0;
            }
            else {
                left_changepoint_position = m_tau[right_cp_index - 1].get_position();
            }
            // i is the first bin association bin point that lies within the range of this interval
            unsigned long int i, i_limit;
            if ((left_changepoint_position * number_of_association_matrix_bins) % m_end == 0) {
                i = (left_changepoint_position * number_of_association_matrix_bins) / m_end;
            }
            else {
                i = (left_changepoint_position * number_of_association_matrix_bins) / m_end + 1;
            }
            // now feed in valid values of i to regime_bin_indices
            if ((right_changepoint_position * number_of_association_matrix_bins) % m_end == 0) {
                i_limit = (right_changepoint_position * number_of_association_matrix_bins) / m_end;
            }
            else {
                i_limit = (right_changepoint_position * number_of_association_matrix_bins) / m_end + 1;
            }
            while (i < i_limit) {
                regime_bin_indices.push_back(i);
                i++;
            }
        }
        size_t number_of_bin_indices = regime_bin_indices.size();
        for (unsigned int bin_index = 0; bin_index < number_of_bin_indices; bin_index++) {
            for (unsigned int bin_index_1 = 0; bin_index_1 < number_of_bin_indices; bin_index_1++) {
                unsigned long int i_0 = regime_bin_indices[bin_index], i_1 = regime_bin_indices[bin_index_1];
                association_matrix[i_0][i_1] += 1;
            }
        }
    }
}

// returns the log of the full I prior ratio if we were adding the changepoint that we are proposing removing with regime proposed_regime. Index gives the position of the changepoint that we are proposing to remove. new_regime gives whether the proposed_regime is a new regime, either in the sense that if the changepoint were deleted it would remove regime proposed_regime, or in the sense that if we were proposing a changepoint here we would be proposing a new regime for it. actual_regime gives the actual regime of the changepoint that we are removing
double particle::full_log_I_prior_smc_remove_ratio(const int & index, const unsigned int & process, const unsigned int & previous_regime, const unsigned int & proposed_regime, const bool & new_regime, const unsigned int & actual_regime, const bool & removing_unobserved_regime, const unsigned int & subsequent_regime) {
    double q;
    // check if the number of regimes has been reduced before running this, if so then reduce the number of regimes by 1.
    double number_of_regimes = static_cast< double >(get_number_of_regimes(process)) - (removing_unobserved_regime ? 1.0 : 0.0);
    if (subsequent_regime == numeric_limits< unsigned int >::max()) {
        if (new_regime) {
            q = log(m_dirichlet_alpha);
            double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions() - 1);
            q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime) - gsl_sf_lngamma(m_dirichlet_alpha + 1 + number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
            q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
            for (unsigned int regime = 0; regime < previous_regime; regime++) {
                if (!removing_unobserved_regime || (regime != actual_regime)) { // only do this if we haven't removed this regime - if actual_regime is removed, then actual_regime can't be used to compute q.
                    double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                    q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                    q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
                }
            }
            for (unsigned int regime = previous_regime + 1; regime < get_number_of_regimes(process); regime++) {
                if (!removing_unobserved_regime || (regime != actual_regime)) { // only do this if we haven't removed this regime
                    double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                    q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                    q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
                }
            }
            q += log(1 - m_rho);
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
        }
        else {
            double transitions_from_previous_to_proposed = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(proposed_regime)) - (proposed_regime == actual_regime ? 1 : 0);
            double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions() - 1);
            q = log(m_dirichlet_alpha + transitions_from_previous_to_proposed) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
        }
    }
    else {
        if (new_regime) {
            double transitions_from_previous_to_subsequent_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(subsequent_regime));
            if (actual_regime == subsequent_regime) {
                if (previous_regime == actual_regime) {
                    transitions_from_previous_to_subsequent_regime--;
                }
            }
            else {
                if (previous_regime != actual_regime) {
                    transitions_from_previous_to_subsequent_regime++;
                }
            }
            q = 2 * log(m_dirichlet_alpha) - log(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - log(m_dirichlet_alpha + transitions_from_previous_to_subsequent_regime - 1);
            for (unsigned int regime = 0; regime < get_number_of_regimes(process); regime++) {
                if (!removing_unobserved_regime || (regime != actual_regime)) { // only do this if we haven't removed this regime - if actual_regime is removed, then actual_regime can't be used to compute q.
                    double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                    transitions_out_of_regime -= (actual_regime == regime ? 1 : 0);
                    q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) + gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
                }
            }
            q += log(1 - m_rho);
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
        }
        else {
            if (proposed_regime == previous_regime) {
                double transitions_from_previous_to_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(previous_regime));
                double transitions_out_of_previous_regime = static_cast< double >(m_regimes[process][previous_regime].get_number_of_right_transitions());
                if (previous_regime == actual_regime) {
                    transitions_from_previous_to_previous_regime--;
                    transitions_out_of_previous_regime--;
                }
                else {
                    if (previous_regime == subsequent_regime) {
                        transitions_from_previous_to_previous_regime++;
                    }
                }
                q = log(m_dirichlet_alpha + transitions_from_previous_to_previous_regime) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_previous_regime);
            }
            else {
                if (subsequent_regime == proposed_regime) {
                    double transitions_from_subsequent_to_subsequent_regime = static_cast< double >(m_regimes[process][subsequent_regime].get_right_transitions_histogram_element(subsequent_regime));
                    double transitions_out_of_subsequent_regime = static_cast< double >(m_regimes[process][subsequent_regime].get_number_of_right_transitions());
                    if (actual_regime == subsequent_regime) {
                        transitions_out_of_subsequent_regime--;
                        transitions_from_subsequent_to_subsequent_regime--;
                    }
                    else {
                        if (previous_regime == subsequent_regime) {
                            transitions_from_subsequent_to_subsequent_regime++;
                        }
                    }
                    q = log(m_dirichlet_alpha + transitions_from_subsequent_to_subsequent_regime) - log(m_dirichlet_alpha * number_of_regimes + transitions_out_of_subsequent_regime);
                }
                else { // know proposed != previous and proposed != subsequent
                    // therefore these three transition counts are definitely different
                    double transitions_from_previous_to_proposed_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(proposed_regime));
                    double transitions_from_proposed_to_subsequent_regime = static_cast< double >(m_regimes[process][proposed_regime].get_right_transitions_histogram_element(subsequent_regime));
                    double transitions_from_previous_to_subsequent_regime = static_cast< double >(m_regimes[process][previous_regime].get_right_transitions_histogram_element(subsequent_regime));
                    if (actual_regime == subsequent_regime) {
                        if (previous_regime == actual_regime) {
                            transitions_from_previous_to_subsequent_regime--;
                        }
                        else {
                            if (proposed_regime == actual_regime) {
                                transitions_from_proposed_to_subsequent_regime--;
                            }
                        }
                    }
                    else {
                        if (previous_regime != actual_regime) {
                            if (proposed_regime == actual_regime) {
                                transitions_from_previous_to_proposed_regime--;
                                transitions_from_proposed_to_subsequent_regime--;
                            }
                            transitions_from_previous_to_subsequent_regime++;
                        }
                    }
                    double transitions_out_of_proposed_regime = static_cast< double >(m_regimes[process][proposed_regime].get_number_of_right_transitions());
                    if (proposed_regime == actual_regime) {
                        transitions_out_of_proposed_regime--;
                    }
                    q = log(m_dirichlet_alpha + transitions_from_previous_to_proposed_regime) + log(m_dirichlet_alpha + transitions_from_proposed_to_subsequent_regime) - log(m_dirichlet_alpha + transitions_from_previous_to_subsequent_regime - 1) - log(number_of_regimes * m_dirichlet_alpha + transitions_out_of_proposed_regime);
                }
            }
        }
    }
    return q;
}

// following_regime and proposed_regime default to -1. Gives the log of the full I prior ratio for adding the separator as a changepoint after it has been removed.
double particle::log_smc_resampling_separator_changepoint_prior_ratio(const unsigned int & process, const bool & no_following_regime, const bool & new_regime, const bool & removing_unobserved_regime, const unsigned int & actual_regime, const unsigned int & following_regime, const unsigned int & proposed_regime) {
    double q = 0;
    double number_of_regimes = static_cast< double >(get_number_of_regimes(process)) - (removing_unobserved_regime ? 1.0 : 0.0);
    if (no_following_regime) {
        if (new_regime) {
            for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
                if (!removing_unobserved_regime || (actual_regime != regime)) {
                    double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions());
                    q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                    q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
                }
            }
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
            q += log(1 - m_rho);
        }
    }
    else {
        if (new_regime) {
            for (unsigned int regime = 0; regime < m_regimes[process].size(); regime++) {
                if (!removing_unobserved_regime || (actual_regime != regime)) {
                    double transitions_out_of_regime = static_cast< double >(m_regimes[process][regime].get_number_of_right_transitions() - (actual_regime == regime ? 1 : 0));
                    q += gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime) - gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha + transitions_out_of_regime);
                    q += gsl_sf_lngamma(m_dirichlet_alpha + number_of_regimes * m_dirichlet_alpha) - gsl_sf_lngamma(number_of_regimes * m_dirichlet_alpha);
                }
            }
            q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
            q += log(m_dirichlet_alpha) - log(number_of_regimes * m_dirichlet_alpha + m_dirichlet_alpha);
            q += log(1 - m_rho);
        }
        else {
            if (actual_regime == proposed_regime) {
                double transitions_from_proposed_to_following_regime = static_cast< double >(m_regimes[process][proposed_regime].get_right_transitions_histogram_element(following_regime));
                double transitions_out_of_proposed_regime = static_cast< double >(m_regimes[process][proposed_regime].get_number_of_right_transitions());
                q += log(transitions_from_proposed_to_following_regime - 1 + m_dirichlet_alpha) - log(transitions_out_of_proposed_regime - 1 + (number_of_regimes * m_dirichlet_alpha));
            }
            else {
                double transitions_from_proposed_to_following_regime = static_cast< double >(m_regimes[process][proposed_regime].get_right_transitions_histogram_element(following_regime));
                double transitions_out_of_proposed_regime = static_cast< double >(m_regimes[process][proposed_regime].get_number_of_right_transitions());
                q += log(transitions_from_proposed_to_following_regime + m_dirichlet_alpha) - log(transitions_out_of_proposed_regime + number_of_regimes * m_dirichlet_alpha);
            }
        }
    }
    return q;
}

void particle::remove_full_changepoint_smc(const unsigned int & index, const vector< vector< double > > & right_sufficient_statistics_rev, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_rev, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_rev, const vector< double > & number_of_observations_rev, const vector< bool > & removing_unobserved_regimes) {
	// get the regime indices for the index'th changepoint
	vector< unsigned int > removed_regimes(0);
	removed_regimes = m_tau[index].get_full_index();

	// remove the changepoint from m_tau
	m_tau.erase(m_tau.begin() + index);
	m_dimension--;
	size_t number_of_processes = removed_regimes.size();
	if (m_include_separator) {
		// if there are separators, decrease the index of any separator greater than index
		decrease_separator_indices_greater_than_index(index);
	}

	for (unsigned int j = 0; j < number_of_processes; j++) {
		unsigned int prev_regime_index;
		// calculate the regime index of the previous changepoint
		if (index == 0) {
			prev_regime_index = 0;
		}
		else {
			prev_regime_index = m_tau[index - 1].get_full_index_row(j);
		}

		// remove the interval and transition out of the removed_regime
		// replace the right transition from previous regime to removed_regime with a transition to subsequent_regime (unless new_regime is equal to subsequent_regime), removing a transition out of previous regime if the changepoint is removed from the end of m_tau.
		// calculate m_trace_index for the changepoint we are removing
		if (is_changepoint_index_separator_index(index) || index == m_dimension) {
			if (index == m_dimension) {
				m_trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
			}
			m_regimes[j][prev_regime_index].remove_right_transition(index);
			m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
		}
		else {
			unsigned int subs_regime_index = m_tau[index].get_full_index_row(j);
			if (removed_regimes[j] != subs_regime_index) {
				m_regimes[j][prev_regime_index].alter_right_transition(index, subs_regime_index);
			}
			m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
		}
		// calculate the trace index
		is_changepoint_index_separator_index(index - 1); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index) has just been run and index now corresponds to the next changepoint

		// remove right sufficient statistics from the regime for the deleted cp and add them to the previous regime (unless both regimes are equal). Calculate the log likelihood for each regime
		if (removed_regimes[j] != prev_regime_index) {
			m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_rev[j]);
			m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_log_likelihoods_without_right_sufficient_statistics_rev[j]);
			m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, number_of_observations_rev[j]);
			m_regimes[j][prev_regime_index].add_sufficient_statistics(right_sufficient_statistics_rev[j]);
			m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_right_sufficient_statistics_rev[j]);
			m_regimes[j][prev_regime_index].add_observations(m_trace_index, number_of_observations_rev[j]);
		}

		if (removing_unobserved_regimes[j]) {
			// if we are removing the top regime, then this regime must be deleted
			if (!m_unobserved_regimes[j][removed_regimes[j]]) {
				cerr << "regime to delete is not unobserved" << endl;
			}
			m_regimes[j].erase(m_regimes[j].begin() + removed_regimes[j]);
			// for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
			for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
				m_regimes[j][regime].remove_unobserved_regime(removed_regimes[j]);
			}
            for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
                m_tau[cp_index].decrease_full_index_row_for_removing_unobserved_regime_if_necessary(j, removed_regimes[j]);
            }
			m_unobserved_regimes[j].erase(m_unobserved_regimes[j].begin() + removed_regimes[j]);
			m_number_of_unobserved_regimes[j]--;
		}

		// decrease the values for right indices in every regime so that index + 2 -> index + 1, index + 3 -> index + 2, ...
		for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
			m_regimes[j][regime].decrease_right_indices_greater_than(index + 1);
		}
	}
}

void particle::move_full_changepoint_smc(const unsigned int & index, const changepoint & new_changepoint, const vector< unsigned int > & new_regimes, const bool & tau_h_greater_than_tau_h_prime, const vector< vector< double > > & right_sufficient_statistics, const vector< double > & right_number_of_observations, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & middle_sufficient_statistics, const vector< double > & middle_number_of_observations, const vector< double > & previous_log_likelihoods_without_right_sufficient_statistics, const vector< double > & previous_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & previous_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & previous_log_likelihoods_without_middle_sufficient_statistics, const vector< vector< double > > & log_likelihoods_with_right_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_with_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_middle_sufficient_statistics, const vector< double > & actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes) {
	vector< unsigned int > removed_regimes(0);
	removed_regimes = m_tau[index].get_full_index();

	m_tau[index] = new_changepoint;
	m_tau[index].set_full_index(new_regimes);
	size_t number_of_processes = removed_regimes.size();
	for (unsigned int j = 0; j < number_of_processes; j++) {
		unsigned int prev_regime_index;
		// calculate the regime index of the previous changepoint
		if (index == 0) {
			prev_regime_index = 0;
		}
		else {
			prev_regime_index = m_tau[index - 1].get_full_index_row(j);
		}

		// if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
		bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
		if (adding_new_regime) {
			for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
				// increase the size of the transitions histogram for each of the regimes in this process.
				m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
			}
			m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(right_sufficient_statistics[j].size(), 0), m_number_of_traces, m_trace_index, 0));
			// the interval will be added later
			m_unobserved_regimes[j].push_back(true);
			m_number_of_unobserved_regimes[j]++;
		}

		// remove right sufficient statistics from the regime for the deleted cp and add them to the previous regime (unless both regimes are equal). Calculate the log likelihood for each regime
		if (removed_regimes[j] == prev_regime_index) {
			if (prev_regime_index != new_regimes[j]) {
				m_regimes[j][prev_regime_index].remove_sufficient_statistics(right_sufficient_statistics[j]);
				m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_without_right_sufficient_statistics[j]);
				m_regimes[j][prev_regime_index].remove_observations(m_trace_index, right_number_of_observations[j]);
				m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics[j]);
				m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]]);
				m_regimes[j][new_regimes[j]].add_observations(m_trace_index, right_number_of_observations[j]);
			}
		}
		else {
			if (new_regimes[j] == prev_regime_index) {
				m_regimes[j][prev_regime_index].add_sufficient_statistics(right_sufficient_statistics_reverse[j]);
				m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_right_sufficient_statistics_reverse[j]);
				m_regimes[j][prev_regime_index].add_observations(m_trace_index, right_number_of_observations_reverse[j]);
				m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
				m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
				m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, right_number_of_observations_reverse[j]);
			}
			else {
				if (new_regimes[j] == removed_regimes[j]) {
					if (tau_h_greater_than_tau_h_prime) {
						m_regimes[j][prev_regime_index].remove_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_without_middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].remove_observations(m_trace_index, middle_number_of_observations[j]);
						m_regimes[j][removed_regimes[j]].add_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_with_middle_sufficient_statistics[j]);
						m_regimes[j][removed_regimes[j]].add_observations(m_trace_index, middle_number_of_observations[j]);
					}
					else {
						m_regimes[j][prev_regime_index].add_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].add_observations(m_trace_index, middle_number_of_observations[j]);
						m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_middle_sufficient_statistics[j]);
						m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, middle_number_of_observations[j]);
					}
				}
				else {
					if (tau_h_greater_than_tau_h_prime) {
						m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
						m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, right_number_of_observations_reverse[j]);
						m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics[j]);
						m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]]);
						m_regimes[j][new_regimes[j]].add_observations(m_trace_index, right_number_of_observations[j]);
						m_regimes[j][prev_regime_index].remove_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_without_middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].remove_observations(m_trace_index, middle_number_of_observations[j]);
					}
					else {
						m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
						m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_regime_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
						m_regimes[j][removed_regimes[j]].remove_observations(m_trace_index, right_number_of_observations_reverse[j]);
						m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics[j]);
						m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods_with_right_sufficient_statistics[j][new_regimes[j]]);
						m_regimes[j][new_regimes[j]].add_observations(m_trace_index, right_number_of_observations[j]);
						m_regimes[j][prev_regime_index].add_sufficient_statistics(middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].set_log_likelihood(previous_log_likelihoods_with_middle_sufficient_statistics[j]);
						m_regimes[j][prev_regime_index].add_observations(m_trace_index, middle_number_of_observations[j]);
					}
				}
			}
		}

		// replace the interval for removed_regime with an interval for new_regime (if they differ)
		if (index == m_dimension - 1 || is_changepoint_index_separator_index(index + 1)) {
			if (new_regimes[j] != removed_regimes[j]) {
				m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
				m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);
				m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
			}
		}
		else {
			unsigned int subs_regime_index = m_tau[index + 1].get_full_index_row(j);
			if (removed_regimes[j] != new_regimes[j]) {
				m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
				m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j], subs_regime_index);
				m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
			}
		}
		// calculate the trace index
		is_changepoint_index_separator_index(index); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index + 1) has just been run

		// if we are removing an unobserved regime and we don't add it back in, then the top regime must be removed
		if (removing_unobserved_regimes[j]) {
			if (!m_unobserved_regimes[j][removed_regimes[j]]) {
				cerr << "regime to delete is not unobserved" << endl;
			}
			m_regimes[j].erase(m_regimes[j].begin() + removed_regimes[j]);
			// for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
			for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
				m_regimes[j][regime].remove_unobserved_regime(removed_regimes[j]);
			}
            for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
                m_tau[cp_index].decrease_full_index_row_for_removing_unobserved_regime_if_necessary(j, removed_regimes[j]);
            }
			m_unobserved_regimes[j].erase(m_unobserved_regimes[j].begin() + removed_regimes[j]);
			m_number_of_unobserved_regimes[j]--;
		}
	}
}

void particle::resample_full_changepoint_smc(const unsigned int & index, const vector< unsigned int > & new_regimes, const vector< vector< double > > & right_sufficient_statistics_reverse, const vector< double > & right_number_of_observations_reverse, const vector< vector< double > > & regime_log_likelihoods_with_right_sufficient_statistics_reverse, const vector< double > & actual_log_likelihoods_without_right_sufficient_statistics_reverse, const vector< bool > & removing_unobserved_regimes) {
	// calculate the regimes that we are going to lose, and set m_tau[index] to be associated with the new regimes.
	vector< unsigned int > removed_regimes(0);
	removed_regimes = m_tau[index].get_full_index();
	m_tau[index].set_full_index(new_regimes);
	size_t number_of_processes = new_regimes.size();

	for (unsigned int j = 0; j < number_of_processes; j++) {
		// if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
		bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
		if (adding_new_regime) {
			for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
				// increase the size of the transitions histogram for each of the regimes in this process.
				m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
			}
			m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(right_sufficient_statistics_reverse[j].size(), 0), m_number_of_traces, m_trace_index, 0));
			m_unobserved_regimes[j].push_back(true);
			m_number_of_unobserved_regimes[j]++;
		}

		unsigned int prev_regime_index;
		// calculate the regime index of the previous changepoint
		if (index == 0) {
			prev_regime_index = 0;
		}
		else {
			prev_regime_index = m_tau[index - 1].get_full_index_row(j);
		}

		// replace the interval for removed_regime with an interval for new_regime (if they differ)
		if (index == m_dimension - 1 || is_changepoint_index_separator_index(index + 1)) {
			if (index == m_dimension - 1) {
				// there is a separator changepoint at the beginning of the last trace, so we must be in the last trace.
				m_h_trace_index = static_cast< unsigned int >(m_number_of_traces - 1);
			}
			else {
				// calculate the trace index
				is_changepoint_index_separator_index(index); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index + 1) has just been run
				m_h_trace_index = m_trace_index;
			}
			if (new_regimes[j] != removed_regimes[j]) {
				if (!is_changepoint_index_separator_index(index)) {
					m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
				}
				m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);
				m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
			}
		}
		else {
			m_h_trace_index = m_trace_index;
			unsigned int subs_regime_index = m_tau[index + 1].get_full_index_row(j);
			if (removed_regimes[j] != new_regimes[j]) {
				m_regimes[j][removed_regimes[j]].remove_interval(index + 1, m_unobserved_regimes[j], removed_regimes[j], m_number_of_unobserved_regimes[j]);
				m_regimes[j][new_regimes[j]].add_interval(index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j], subs_regime_index);
				if (!is_changepoint_index_separator_index(index)) {
					m_regimes[j][prev_regime_index].alter_right_transition(index, new_regimes[j]);
				}
			}
		}

		if (removed_regimes[j] != new_regimes[j]) {
			m_regimes[j][removed_regimes[j]].remove_sufficient_statistics(right_sufficient_statistics_reverse[j]);
			m_regimes[j][removed_regimes[j]].set_log_likelihood(actual_log_likelihoods_without_right_sufficient_statistics_reverse[j]);
			m_regimes[j][removed_regimes[j]].remove_observations(m_h_trace_index, right_number_of_observations_reverse[j]);
			m_regimes[j][new_regimes[j]].add_sufficient_statistics(right_sufficient_statistics_reverse[j]);
			m_regimes[j][new_regimes[j]].set_log_likelihood(regime_log_likelihoods_with_right_sufficient_statistics_reverse[j][new_regimes[j]]);
			m_regimes[j][new_regimes[j]].add_observations(m_h_trace_index, right_number_of_observations_reverse[j]);
		}

		// if we are removing an unobserved regime it must be removed
		if (removing_unobserved_regimes[j]) {
			if (!m_unobserved_regimes[j][removed_regimes[j]]) {
				cerr << "regime to delete is not unobserved" << endl;
			}
			m_regimes[j].erase(m_regimes[j].begin() + removed_regimes[j]);
			// for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
			for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
				m_regimes[j][regime].remove_unobserved_regime(removed_regimes[j]);
			}
            for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
                m_tau[cp_index].decrease_full_index_row_for_removing_unobserved_regime_if_necessary(j, removed_regimes[j]);
            }
			m_unobserved_regimes[j].erase(m_unobserved_regimes[j].begin() + removed_regimes[j]);
			m_number_of_unobserved_regimes[j]--;
		}
	}
}

double particle::full_log_I_prior_add_ratio_always_new_regime(const unsigned int & process, const bool & new_regime) {
    double number_of_regimes = static_cast< double >(get_number_of_regimes(process));
    if (new_regime) {
        double q = log(1 - m_rho);
        q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
        return q;
    }
    else {
        return 0;
    }
}

double particle::full_log_I_prior_remove_ratio_always_new_regime(const unsigned int & process, const bool & new_regime, const bool & removing_unobserved_regime) {
    double q;
    // check if the number of regimes has been reduced before running this, if so then reduce the number of regimes by 1.
    double number_of_regimes = static_cast< double >(get_number_of_regimes(process)) - (removing_unobserved_regime ? 1.0 : 0.0);
    if (new_regime) {
        q = log(1 - m_rho);
        q += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
    }
    else {
        q = 0;
    }
    return q;
}

double particle::calculate_and_get_add_unobserved_regimes_full_I_prior_ratio_always_new_regime(const unsigned int & process) {
    double ratio = 0;
    double number_of_regimes = static_cast< double >(m_regimes[process].size());
    ratio += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
    return ratio;
}

double particle::calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio_always_new_regime(const unsigned int & process) {
    double ratio = 0;
    double number_of_regimes = static_cast< double >(m_regimes[process].size());
    ratio += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes - 1));
    return ratio;
}

#endif

/*void particle::add_new_interval_changepoints(const vector< changepoint > & new_changepoints, const vector< vector< vector< double > > > & sufficient_statistics, const vector< vector< double > > & log_likelihoods, const vector< vector< double > > & number_of_observations) { // sufficient_statistics has a vector of vectors of doubles for each changepoint
size_t number_of_processes = sufficient_statistics[0].size();
m_trace_index = static_cast< unsigned int >(m_separators.size());
for (unsigned int cp_count = 0; cp_count < new_changepoints.size(); cp_count++) {
unsigned int new_cp_index = m_dimension;
// insert new_changepoint into m_tau
m_tau.insert(m_tau.begin() + new_cp_index, new_changepoints[cp_count]);
m_dimension++;
m_log_k_prior += calculate_and_get_add_cp_k_prior_ratio();
// new_regimes is a vector of the regimes that will be adopted at the spearator changepoint
vector< unsigned int > new_regimes = vector< unsigned int >(number_of_processes);
for (unsigned int j = 0; j < number_of_processes; j++) {
new_regimes[j] = static_cast< unsigned int >(m_regimes[j].size());
}

for (unsigned int j = 0; j < number_of_processes; j++) {
unsigned int prev_regime_index;
// calculate the regime index of the previous changepoint
if (new_cp_index == 0) {
prev_regime_index = 0;
}
else {
prev_regime_index = m_tau[new_cp_index - 1].get_full_index_row(j);
}
// if adding a new regime, insert it into the collection of existing regimes, and change the other regimes so that they know the number of regimes has increased.
bool adding_new_regime = new_regimes[j] == m_regimes[j].size();
if (adding_new_regime) {
for (unsigned int i = 0; i < m_regimes[j].size(); i++) {
// increase the size of the transitions histogram for each of the regimes in this process.
m_regimes[j][i].add_regime_to_transitions_histogram(new_regimes[j]);
}
// insert new regime into m_regimes[j]. We don't include the right changepoint position, right transition, sufficient statistics or number_of_observations here as it will be added later
// these will be converted to the correct values later
m_regimes[j].push_back(regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size() + 1, 0), vector< double >(sufficient_statistics[j].size(), 0), m_number_of_traces, m_trace_index, 0));
// these unobserved regime values will be corrected later
m_unobserved_regimes[j].push_back(true);
m_number_of_unobserved_regimes[j]++;
m_log_regimes_prior += log(1 - m_rho);
}

// increase the values for the right indices in every regime so that index + 1 -> index + 2, ...
for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
m_regimes[j][regime].increase_right_indices_greater_than_or_equal_to(new_cp_index + 1);
}

// insert the interval and transition out of the new_regime
// replace the right transition from previous regime to subsequent regime with a transition to new_regime (unless new_regime is equal to subsequent_regime), adding a transition out of previous regime if the new changepoint is inserted at the end of m_tau.
m_regimes[j][prev_regime_index].alter_right_transition(new_cp_index, new_regimes[j]);
m_regimes[j][new_regimes[j]].add_interval(new_cp_index + 1, m_unobserved_regimes[j], new_regimes[j], m_number_of_unobserved_regimes[j]);

// assign right sufficient statistics and observation numbers to the new regime and take them away from the previous regime (unless both regimes are equal). Calculate the log likelihood for each regime
m_regimes[j][prev_regime_index].remove_sufficient_statistics(sufficient_statistics[cp_count][j]);
m_log_likelihood += log_likelihoods[cp_count][j] - m_regimes[j][prev_regime_index].get_log_likelihood(); // this will be zero except for the first time
m_regimes[j][prev_regime_index].set_log_likelihood(log_likelihoods[cp_count][j]);
m_regimes[j][prev_regime_index].remove_observations(m_trace_index, number_of_observations[cp_count][j]);
m_regimes[j][new_regimes[j]].add_sufficient_statistics(sufficient_statistics[cp_count][j]);
m_regimes[j][new_regimes[j]].set_log_likelihood(log_likelihoods[cp_count + 1][j]);
m_log_likelihood += m_regimes[j][new_regimes[j]].get_log_likelihood();
m_regimes[j][new_regimes[j]].add_observations(m_trace_index, number_of_observations[cp_count][j]);
}
}
m_log_full_I_prior = calculate_and_get_log_full_I_prior(number_of_processes);
}

*/

// finds where a new changepoint would go
/*unsigned int particle::index_finder(const changepoint & new_changepoint, unsigned int lower_bound, unsigned int upper_bound){
 if(lower_bound == m_dimension)
 return m_dimension;
 if(m_dimension == 0 || new_changepoint < m_tau[lower_bound])
 return lower_bound;
 if(!upper_bound)
 upper_bound = m_dimension - 1;
 if(new_changepoint > m_tau[upper_bound])
 return upper_bound + 1;
 
 unsigned int m = (lower_bound + upper_bound) / 2;
 while(upper_bound - lower_bound > 1){
 (new_changepoint < m_tau[m]) ? upper_bound = m : lower_bound = m;
 m = (lower_bound + upper_bound) / 2;
 }
 return m + 1;
 }*/

/*
 changepoint particle::get_binary_left_cp(const size_t & process, const unsigned int & index){
 if (index == 0){
 return m_intercept_changepoint;
 } else {
 int left_index = get_binary_left_index(process, index);
 if (left_index == -1) {
 return m_intercept_changepoint;
 }
 else {
 return m_tau[m_binaries[process][m_tau[index - 1].get_binary_index_row(process)].get_left_index()];
 }
 }
 }
 
 changepoint particle::get_binary_right_cp(const size_t & process, const unsigned int & index){
 if (index == 0){
 return m_tau[m_binaries[process][1].get_left_index()];
 } else {
 return m_tau[m_binaries[process][m_tau[index - 1].get_binary_index_row(process) + 1].get_left_index()];
 }
 }*/

/*void particle::increase_log_full_I_prior_unobserved(const double & log_full_I_prior_ratio, const vector< bool > & unobserved_regime_change, const bool & adding_unobserved_regimes) {
 m_log_full_I_prior += log_full_I_prior_ratio;
 if (adding_unobserved_regimes) {
 for (unsigned int process = 0; process < unobserved_regime_change.size(); process++) {
 if (unobserved_regime_change[process]) {
 double number_of_regimes = static_cast< double >(m_regimes[process].size());
 m_log_full_I_prior -= static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
 m_log_full_separators_prior += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes + 1));
 }
 }
 }
 else {
 for (unsigned int process = 0; process < unobserved_regime_change.size(); process++) {
 if (unobserved_regime_change[process]) {
 double number_of_regimes = static_cast< double >(m_regimes[process].size());
 m_log_full_I_prior -= static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes - 1));
 m_log_full_separators_prior += static_cast< double >(m_separator_indices.size()) * log(number_of_regimes / (number_of_regimes - 1));
 }
 }
 }
 }*/


/*void particle::add_unobserved_regimes(const vector< bool > & adding_unobserved_regimes, const size_t & number_of_processes) {
 for (unsigned int j = 0; j < number_of_processes; j++) {
 if (adding_unobserved_regimes[j]) {
 // choose the new unobserved regime (the only restriction is that it can't be regime 0)
 unsigned int new_unobserved_regime = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_regimes[j].size()) + 1);
 size_t size_of_sufficient_statistics = m_regimes[j][0].get_size_of_sufficient_statistics();
 // insert the new unobserved regime into m_regimes[j]
 m_regimes[j].insert(m_regimes[j].begin() + new_unobserved_regime, regime(vector< int >(0), vector< unsigned int >(0), vector< unsigned int >(m_regimes[j].size(), 0), vector< double >(size_of_sufficient_statistics, 0), m_number_of_traces, 0, 0));
 // set the likelihood to be 0 for this new regime
 m_regimes[j][new_unobserved_regime].set_log_likelihood(0);
 // for all other regimes, include an empty transition in m_right_transitions_histogram and, for each element in m_right_transitions, if the element is greater than or equal to the new_unobserved_regime then add 1 to it
 for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
 m_regimes[j][regime].insert_new_unobserved_regime(new_unobserved_regime);
 }
 for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
 m_tau[cp_index].increase_full_index_row_for_inserting_unobserved_regime_if_necessary(j, new_unobserved_regime);
 }
 m_unobserved_regimes[j].insert(m_unobserved_regimes[j].begin() + new_unobserved_regime, true);
 m_number_of_unobserved_regimes[j]++;
 }
 }
 }
 
 void particle::remove_unobserved_regimes(const vector< bool > & removing_unobserved_regimes, const size_t & number_of_processes) {
 for (unsigned int j = 0; j < number_of_processes; j++) {
 if (removing_unobserved_regimes[j]) {
 // choose which unobserved regime to remove using sub-linear method - keep guessing at unobserved regimes
 unsigned int regime_to_remove = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_regimes[j].size() - 1) + 1);
 bool unobserved = m_unobserved_regimes[j][regime_to_remove];
 while (!unobserved) {
 regime_to_remove = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_regimes[j].size() - 1) + 1);
 unobserved = m_unobserved_regimes[j][regime_to_remove];
 }
 m_regimes[j].erase(m_regimes[j].begin() + regime_to_remove);
 // for all other regimes, remove any trace of this unobserved regime (i.e. delete the 0 from m_right_transitions_histogram, decrease the indices of transitions greate than regime_to_remove by 1
 for (unsigned int regime = 0; regime < m_regimes[j].size(); regime++) {
 m_regimes[j][regime].remove_unobserved_regime(regime_to_remove);
 }
 for (unsigned int cp_index = 0; cp_index < m_dimension; cp_index++) {
 m_tau[cp_index].decrease_full_index_row_for_removing_unobserved_regime_if_necessary(j, regime_to_remove);
 }
 m_unobserved_regimes[j].erase(m_unobserved_regimes[j].begin() + regime_to_remove);
 m_number_of_unobserved_regimes[j]--;
 }
 }
 }*/
