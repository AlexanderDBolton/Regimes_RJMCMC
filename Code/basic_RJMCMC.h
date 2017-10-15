
#ifndef BASIC_RJMCMC_H
#define BASIC_RJMCMC_H

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "basic_particle.h"
#include "basic_changepoint.h"
#include "multiple_processes_regime.h"

class rj{
public:
    rj(const unsigned long int & start_time = 0, const unsigned long int & end_time = 1, const double & p = 1, const double & var_p = 0, const unsigned long int & burnin = 0, const unsigned long int & iterations = 1, const unsigned long int & thinning = 1, const unsigned long int & number_of_association_matrix_bins = 0, const vector< unsigned long int > & starting_changepoints = vector< unsigned long int >(0), const vector< unsigned long int > & separators = vector< unsigned long int >(0), const vector< unsigned long int > & trace_lengths = vector< unsigned long int >(0), const unsigned int & seed = 1, mult_process * pm_ptr = NULL, const unsigned long int & intercept = 0, const unsigned long int & diff = 1); // Constructor for normal RJMCMC
    void set_binary_burnin_iterations_thinning( const unsigned long int & binary_burnin, const unsigned long int & binary_iterations, const unsigned long int & binary_thinning ) {m_binary_burnin = binary_burnin; m_binary_iterations = binary_iterations; m_binary_thinning = binary_thinning; }
    void set_full_burnin_iterations_thinning( const unsigned long int & full_burnin, const unsigned long int & full_iterations, const unsigned long int & full_thinning ) {m_full_burnin = full_burnin; m_full_iterations = full_iterations; m_full_thinning = full_thinning; }
    void record_basic_samples(const bool & record = false, const unsigned int & number_of_bins = 100) { m_recording_basic_samples = record; m_number_of_changepoint_bins = number_of_bins; m_basic_MAP_log_posterior = -1e300; }
    void record_binary_samples(const bool & record = false, const unsigned int & number_of_bins = 100) { m_recording_binary_samples = record; m_number_of_changepoint_bins = number_of_bins; m_binary_MAP_log_posterior = -1e300; }
    void record_full_samples(const bool & record = false, const string & data_file = "", const unsigned int & number_of_bins = 100) { m_recording_full_samples = record; m_data_file = data_file; m_number_of_changepoint_bins = number_of_bins; m_full_MAP_log_posterior = m_particle.get_full_log_posterior(); m_full_MAP_particle = m_particle; m_full_MAP_dimension = m_particle.get_dimension(); m_recorded_full_birth_proposals = 0; m_recorded_full_birth_acceptances = 0; m_recorded_full_death_proposals = 0; m_recorded_full_death_acceptances = 0; m_recorded_full_move_proposals = 0; m_recorded_full_move_acceptances = 0; m_recorded_full_resample_proposals = 0; m_recorded_full_resample_acceptances = 0; m_recorded_full_unobserveds_proposals = 0; m_recorded_full_unobserveds_acceptances = 0;}
    double calculate_total_basic_log_likelihood();
    double calculate_full_log_acceptance_probability(const double & u1);
    double calculate_process_full_log_I_prior(const double & process, const unsigned int & adding_index = -1, const unsigned int & proposed_regime = -1);
    double calculate_process_full_log_I_prior_without_cp(const double & process, const unsigned int & removing_cp, const unsigned int & proposed_regime = -1);
    double calculate_full_log_I_prior(const vector< unsigned int > & regime_vector, const size_t & number_of_regimes);
    vector< vector< double > > generate_regime_transition_matrix(const vector< unsigned int > & regime_vector, const size_t & number_of_regimes);
    double calculate_log_of_sum_Bq_ratio(const vector< double > & new_I_priors, const double & old_I_prior);
    void adding_basic_changepoint_setup(const unsigned int & trace_index);
    void removing_basic_changepoint_setup(const unsigned int & trace_index);
    void moving_basic_changepoint_setup(const unsigned int & trace_index);
	void basic_acceptance_procedure( const double & u1);
	void basic_recording_procedure();
    void run_basic_simulation();
    void write_basic_MAP_changepoints_to_file(const string & basic_MAP_cps_Filename);
    void write_binary_MAP_changepoints_to_file(const string & binary_MAP_cps_Filename);
    void write_full_MAP_changepoints_to_file(const string & full_MAP_cps_Filename);
    void write_basic_dimension_distribution_to_file(const string & basic_dimension_distribution_Filename);
    void write_binary_dimension_distribution_to_file(const string & binary_dimension_distribution_Filename);
    void write_full_dimension_distribution_to_file(const string & full_dimension_distribution_Filename);
    void write_full_effective_dimension_distribution_to_file(const string & full_effective_dimension_distribution_Filename);
    void write_basic_changepoints_distribution_to_file(const string & basic_changepoints_distribution_Filename);
    void write_binary_changepoints_distribution_to_file(const string & binary_changepoints_distribution_Filename);
    void write_full_changepoints_distribution_to_file(const string & full_changepoints_distribution_Filename, const unsigned long int number_of_iterations);
    void write_number_of_regimes_to_file(const string & number_of_regimes_Filename);
    void write_number_of_observed_regimes_to_file(const string & number_of_regimes_Filename);
    void write_basic_log_posterior_trace_to_file(const string & basic_log_posterior_trace_Filename);
    void write_binary_log_posterior_trace_to_file(const string & binary_log_posterior_trace_Filename);
    void write_full_log_posterior_trace_to_file(const string & full_log_posterior_trace_Filename);
    void write_basic_dimension_trace_to_file(const string & dimension_trace_Filename);
    void write_dimension_trace_to_file(const string & dimension_trace_Filename);
    void write_number_of_regimes_trace_to_file(const string & number_of_regimes_trace_Filename);
    void write_full_acceptance_probabilities_to_file(const string & acceptance_probabilities_Filename);
    void write_similarity_matrix_to_file(const string & similarity_matrix_Filename);
    void write_min_proportion_similarity_matrix_to_file(const string & min_proportion_similarity_matrix_Filename);
    void write_similarity_matrices_to_file(const string & similarity_matrices_Filename);
    void write_min_proportion_similarity_matrices_to_file(const string & min_proportion_similarity_matrices_Filename);
    void write_association_matrix_to_file(const string & association_matrix_Filename);
    void convert_basic_particle_to_binary_particle(const double & beta_alpha);
	void set_binary_marked_vectors();
    double calculate_total_binary_log_likelihood();
    void adding_binary_changepoint_setup(const unsigned int & trace_index);
    void removing_binary_changepoint_setup(const unsigned int & trace_index);
    void moving_binary_changepoint_setup(const unsigned int & trace_index);
    void resampling_binary_changepoint_setup(const unsigned int & trace_index);
    void binary_acceptance_procedure(const double & u1);
    void binary_recording_procedure();
    void run_binary_simulation();
    void convert_binary_particle_to_full_particle(const double & dirichlet_alpha, const double & rho);
	void set_full_marked_vectors();// const size_t & number_of_processes, const vector< vector< vector< double > > > & sufficient_statistics, const vector< vector< double > > & number_of_observations);
    void check_total_full_log_likelihood(particle & P);
    void check_adding_changepoint(const unsigned int & add_cp_index);
    void calculate_vector_descending_order(const size_t & length, vector< double > A, vector< unsigned int > & order);
    void adding_full_changepoint_setup(const unsigned int & trace_index);
    void removing_full_changepoint_setup(const unsigned int & trace_index);
    void moving_full_changepoint_setup(const unsigned int & trace_index);
    void resampling_full_changepoint_setup(const unsigned int & trace_index);
    void altering_unobserved_regimes_setup();
    //void remove_unobserved_regimes_setup();
    void full_acceptance_procedure(const double & u1);
    void full_recording_procedure();
    void update_full_MAP();
    void run_full_simulation();
    particle get_particle() const {return m_particle;}
    
protected:
    bool m_recording_basic_samples;
    bool m_recording_binary_samples;
    bool m_recording_full_samples;
    bool m_recording_association_matrix;
    string m_data_file;
	vector< unsigned int > m_recorded_basic_dimensions;
    vector< unsigned int > m_recorded_binary_dimensions;
	vector< unsigned int > m_recorded_full_dimensions;
    vector< unsigned int > m_recorded_full_effective_dimensions;
	vector< unsigned long int > m_recorded_basic_changepoints;
    vector< unsigned long int > m_recorded_binary_changepoints;
    vector< unsigned long int > m_recorded_full_changepoints;
    vector< vector< size_t > > m_recorded_number_of_regimes;
    double m_recorded_full_birth_acceptances;
    double m_recorded_full_birth_proposals;
    double m_recorded_full_death_acceptances;
    double m_recorded_full_death_proposals;
    double m_recorded_full_move_acceptances;
    double m_recorded_full_move_proposals;
    double m_recorded_full_resample_acceptances;
    double m_recorded_full_resample_proposals;
    double m_recorded_full_unobserveds_acceptances;
    double m_recorded_full_unobserveds_proposals;
    vector< vector< unsigned int > > m_recorded_number_of_observed_regimes;
	unsigned long int m_number_of_changepoint_bins;
	vector< double > m_recorded_basic_log_posteriors;
    vector< double > m_recorded_binary_log_posteriors;
    vector< double > m_recorded_full_log_posteriors;
    vector< vector< double > > m_recorded_similarity_matrix;
    vector< vector< double > > m_recorded_min_proportion_similarity_matrix;
    vector< vector< vector< double > > > m_recorded_similarity_matrices;
    vector< vector< vector< double > > > m_recorded_min_proportion_similarity_matrices;
    vector< vector< vector< double > > > m_association_matrices; // for each process, for each pair of points, calculate the number of times that they fall in the same regime.
    vector< unsigned long int > m_observations_in_each_trace;
	particle m_basic_MAP_particle;
    particle m_binary_MAP_particle;
    particle m_full_MAP_particle;
	double m_basic_MAP_log_posterior;
    double m_binary_MAP_log_posterior;
    double m_full_MAP_log_posterior;
	unsigned int m_basic_MAP_dimension;
    unsigned int m_binary_MAP_dimension;
    unsigned int m_full_MAP_dimension;
    unsigned long int m_start;
    unsigned long int m_end;
    double m_p;
    double m_var_p;
    unsigned long int m_basic_burnin;
    unsigned long int m_basic_iterations;
    unsigned long int m_basic_thinning;
    unsigned long int m_binary_burnin;
    unsigned long int m_binary_iterations;
    unsigned long int m_binary_thinning;
    unsigned long int m_full_burnin;
    unsigned long int m_full_iterations;
    unsigned long int m_full_thinning;
    unsigned long int m_number_of_association_matrix_bins;
    vector< unsigned long int > m_separators;
    unsigned long int m_diff;
    size_t m_number_of_traces;
    unsigned int m_seed;
    mult_process * m_pm_ptr;
    size_t m_number_of_processes;
    unsigned long int m_intercept;
    particle m_particle;
    unsigned int m_dimension;
    changepoint m_end_changepoint;
    double m_b_k;
    double m_d_k;
    double m_m_k;
    double m_r_k;
    double m_au_k;
    double m_ru_k;
    unsigned int m_h;
    double m_log_proposal_ratio;
    double m_log_likelihood_ratio;
    double m_log_k_prior_ratio;
    double m_log_full_I_prior_ratio;
	double m_log_regimes_prior_ratio;
    double m_log_acceptance_prob;
    double m_right_log_likelihood;
    double m_left_log_likelihood;
    int m_binary_left_index;
    unsigned int m_binary_right_index;
    int m_full_left_index; //needed? delete if not used
    unsigned int m_full_right_index; //needed? delete if not used
    vector< double > m_binary_right_log_likelihood;
    vector< double > m_binary_left_log_likelihood;
    vector< double > m_binary_merged_log_likelihood;
	vector< double > m_binary_left_log_likelihood_reverse;
	vector< double > m_binary_right_log_likelihood_reverse;
    vector< vector< double > > m_log_B;
	vector< vector< double > > m_log_B_reverse;
    vector< vector< double > > m_log_q;
	vector< vector< double > > m_log_q_reverse;
    vector< vector< double > > m_log_Bq;
	vector< vector< double > > m_log_Bq_reverse;
    vector< vector< unsigned int > > m_log_Bq_descending_order;
    vector< vector< unsigned int > > m_log_Bq_reverse_descending_order;
	vector< vector< double > > m_left_sufficient_statistics; // needed?
    vector< vector< double > > m_left_sufficient_statistics_reverse; // needed?
	vector< vector< double > > m_right_sufficient_statistics;
    vector< vector< double > > m_right_sufficient_statistics_reverse;
	vector< vector< double > > m_middle_sufficient_statistics;
    vector< double > m_previous_log_likelihoods_without_right_sufficient_statistics;
    vector< vector< double > > m_log_likelihoods_with_right_sufficient_statistics;
    vector< double > m_previous_log_likelihoods_with_right_sufficient_statistics_reverse;
    vector< double > m_actual_log_likelihoods_without_right_sufficient_statistics_reverse;
    vector< double > m_previous_log_likelihoods_without_right_sufficient_statistics_reverse;
    vector< vector< double > > m_log_likelihoods_with_right_sufficient_statistics_reverse;
    vector< double > m_previous_log_likelihoods_with_right_sufficient_statistics;
    vector< double > m_actual_log_likelihoods_without_right_sufficient_statistics;
    vector< double > m_previous_log_likelihoods_without_middle_sufficient_statistics;
    vector< double > m_previous_log_likelihoods_with_middle_sufficient_statistics;
    vector< double > m_actual_log_likelihoods_with_middle_sufficient_statistics;
    vector< double > m_actual_log_likelihoods_without_middle_sufficient_statistics;
    vector< int > m_altering_unobserved_regimes;
    vector< bool > m_removing_unobserved_regimes;
    vector< double > m_log_of_sum_Bq;
    vector< double > m_log_of_sum_Bq_reverse;
    double m_merged_log_likelihood;
    unsigned long int m_new_changepoint_position;
	bool m_tau_h_greater_than_tau_h_prime;
    changepoint m_adding_changepoint;
    const gsl_rng_type * r_type;
    gsl_rng * r;
};

rj::rj(const unsigned long int & start_time, const unsigned long int & end_time, const double & p, const double & var_p, const unsigned long int & basic_burnin, const unsigned long int & basic_iterations, const unsigned long int & basic_thinning, const unsigned long int & number_of_association_matrix_bins, const vector< unsigned long int > & starting_changepoints, const vector< unsigned long int > & separators, const vector< unsigned long int > & trace_lengths, const unsigned int & seed, mult_process * pm_ptr, const unsigned long int & intercept, const unsigned long int & diff):m_start(start_time), m_end(end_time), m_p(p), m_var_p(var_p), m_basic_burnin(basic_burnin), m_basic_iterations(basic_iterations), m_basic_thinning(basic_thinning), m_recording_association_matrix(number_of_association_matrix_bins > 0), m_number_of_association_matrix_bins(number_of_association_matrix_bins), m_separators(separators), m_number_of_traces(separators.size() + 1), m_seed(seed), m_pm_ptr(pm_ptr), m_intercept(intercept), m_diff(diff) {
    changepoint intercept_changepoint(m_intercept);
    intercept_changepoint.set_log_likelihood(m_pm_ptr->calculate_log_likelihood(m_intercept, m_end + 1));
    
    // calculate the number of observations in each trace
    if (trace_lengths.size() > 0) {
        m_observations_in_each_trace = trace_lengths;
        if (trace_lengths.size() != separators.size() + 1) {
            cerr << "sizes of the number of separators and the trace lengths don't match" << endl;
        }
    }
    else {
        m_observations_in_each_trace = vector< unsigned long int >(separators.size() + 1, 0);
        if (separators.size() > 0) {
            cout << "Assuming observations at 1, 2, 3, ..., s_1 + 1, s_1 + 2, ... in traces" << endl;
            m_observations_in_each_trace[0] = (separators[0] - 1 - m_start) * m_diff; // as first observation at time 1 and separators[0] gives the time of the last observation.
            for (unsigned int trace_idx = 1; trace_idx < separators.size(); trace_idx++) {
                m_observations_in_each_trace[trace_idx] = (separators[trace_idx] - separators[trace_idx - 1] - 1) * m_diff;
            }
            m_observations_in_each_trace[separators.size()] = (end_time - separators[separators.size() - 1]) * m_diff;
        }
        else {
            m_observations_in_each_trace[0] = (end_time - m_start) * m_diff;
        }
    }
    
    m_number_of_processes = m_pm_ptr->get_number_of_processes();
    // create the association matrix for each process if required
    if (m_recording_association_matrix) {
        m_association_matrices = vector< vector< vector< double > > >(m_number_of_processes, vector< vector< double > >(m_number_of_association_matrix_bins, vector< double >(m_number_of_association_matrix_bins, 0)));
    }
    
    // create the particle
    m_particle = particle(m_start, m_end, m_separators, intercept_changepoint, m_p, m_var_p);
    // add separator changepoints
    if (m_number_of_traces > 1) {
        m_particle.get_changepoint(-1).set_log_likelihood(m_pm_ptr->calculate_log_likelihood(m_intercept, m_separators[0]));
        changepoint new_separator_changepoint(m_separators[0]);
        if (m_number_of_traces == 2) {
            new_separator_changepoint.set_log_likelihood(m_pm_ptr->calculate_log_likelihood(m_separators[0], m_end + 1));
            m_particle.add_separator_changepoint(new_separator_changepoint, 0);
        }
        else {
            new_separator_changepoint.set_log_likelihood(m_pm_ptr->calculate_log_likelihood(m_separators[0], m_separators[1]));
            m_particle.add_separator_changepoint(new_separator_changepoint, 0);
            for (unsigned int separator_index = 1; separator_index < m_separators.size() - 1; separator_index++) {
                changepoint new_separator_changepoint(m_separators[separator_index]);
                new_separator_changepoint.set_log_likelihood(m_pm_ptr->calculate_log_likelihood(m_separators[separator_index], m_separators[separator_index + 1]));
                m_particle.add_separator_changepoint(new_separator_changepoint, separator_index);
            }
            changepoint new_separator_changepoint(m_separators[m_separators.size() - 1]);
            new_separator_changepoint.set_log_likelihood(m_pm_ptr->calculate_log_likelihood(m_separators[m_separators.size() - 1], m_end + 1));
            m_particle.add_separator_changepoint(new_separator_changepoint, static_cast< unsigned int >(m_separators.size() - 1));
        }
    }
    
    // add starting changepoints
    if (0 < starting_changepoints.size()) {
        for (unsigned int ncp = 0; ncp < starting_changepoints.size(); ncp++) {
            if (m_particle.does_changepoint_exist_in_particle(starting_changepoints[ncp])) {
                cerr << "changepoint already exists in the particle" << endl;
            }
            m_new_changepoint_position = starting_changepoints[ncp];
            changepoint new_changepoint(m_new_changepoint_position);
            unsigned int add_cp_index = m_particle.get_add_changepoint_index();
            m_left_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(add_cp_index - 1).get_position(), m_new_changepoint_position);
            if (add_cp_index == m_particle.get_dimension()) {
                m_right_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_new_changepoint_position, m_end + 1);
            } else {
                m_right_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_new_changepoint_position, m_particle.get_changepoint(add_cp_index).get_position());
            }
            m_log_likelihood_ratio = m_left_log_likelihood + m_right_log_likelihood - m_particle.get_changepoint(add_cp_index - 1).get_log_likelihood();
            m_log_k_prior_ratio = m_particle.calculate_and_get_add_cp_k_prior_ratio();
            m_particle.increase_log_likelihood(m_log_likelihood_ratio);
            m_particle.increase_log_k_prior(m_log_k_prior_ratio);
            new_changepoint.set_log_likelihood(m_right_log_likelihood);
            m_particle.add_changepoint(m_particle.get_add_changepoint_index(), new_changepoint);
            m_particle.get_changepoint(m_particle.get_add_changepoint_index() - 1).set_log_likelihood(m_left_log_likelihood);
        }
    }
    m_end_changepoint = changepoint(m_end + 1);
    m_particle.set_log_likelihood(calculate_total_basic_log_likelihood());
    m_particle.calculate_and_set_log_k_prior();
    gsl_rng_env_setup();
    r_type = gsl_rng_default;
    r = gsl_rng_alloc(r_type);
    gsl_rng_set(r, seed);
}

double rj::calculate_total_basic_log_likelihood() {
    if (m_particle.get_dimension() == 0) {
        return m_pm_ptr->calculate_log_likelihood(m_intercept, m_end + 1);
    }
    if (0.001 < abs(m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(-1).get_position(), m_particle.get_changepoint(0).get_position()) - m_particle.get_changepoint(-1).get_log_likelihood())) {
        cerr << "actual - recorded = " << m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(-1).get_position(), m_particle.get_changepoint(0).get_position()) - m_particle.get_changepoint(-1).get_log_likelihood() << endl;
    }
    double total_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(-1).get_position(), m_particle.get_changepoint(0).get_position());
    double cp_index_log_likelihood;
    for (unsigned int cp_index = 0; cp_index < m_particle.get_dimension() - 1; cp_index++) {
        cp_index_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(cp_index).get_position(), m_particle.get_changepoint(cp_index + 1).get_position());
        if (0.001 < abs(cp_index_log_likelihood - m_particle.get_changepoint(cp_index).get_log_likelihood())) {
            cerr << "actual - recorded = " << cp_index_log_likelihood - m_particle.get_changepoint(cp_index).get_log_likelihood() << endl;
        }
        total_log_likelihood += cp_index_log_likelihood;
    }
    total_log_likelihood += m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(m_particle.get_dimension() - 1).get_position(), m_end + 1);
    if (0.001 < abs(m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(m_particle.get_dimension() - 1).get_position(), m_end + 1) - m_particle.get_changepoint(m_particle.get_dimension() - 1).get_log_likelihood())) {
        cerr << "actual - recorded = " << m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(m_particle.get_dimension() - 1).get_position(), m_end + 1) - m_particle.get_changepoint(m_particle.get_dimension() - 1).get_log_likelihood() << endl;
    }
    return total_log_likelihood;
}

// used to check that calculations are correct
// calculate the log acceptance probability of the MH move that is being proposed (can be greater than 0, not doing min(1, x) bit)
double rj::calculate_full_log_acceptance_probability(const double & u1) {
    double full_log_acceptance_probability = 0;
    if (u1 < m_b_k) {
        full_log_acceptance_probability += log(m_p) - log(1 - m_p);
        full_log_acceptance_probability += log(static_cast< double >(m_end - m_dimension)) - log(static_cast< double >(m_dimension + 1));
        full_log_acceptance_probability += m_dimension == 0 ? -log(4.0) : 0.0;
        double old_I_prior_total = 0;
        for (unsigned int i = 0; i < m_number_of_processes; i++) {
            double old_I_prior = calculate_process_full_log_I_prior(i);
            old_I_prior_total += old_I_prior;
            size_t number_of_regimes = m_particle.get_number_of_regimes(i);
            vector< double > new_I_priors = vector< double >(number_of_regimes + 1);
            for (unsigned int reg = 0; reg < number_of_regimes; reg++) {
                new_I_priors[reg] = calculate_process_full_log_I_prior(i, m_particle.get_add_changepoint_index(), reg);
                if ((reg == number_of_regimes - 1) && (m_particle.is_regime_unobserved(i, number_of_regimes - 1))) {
                    new_I_priors[reg] = -1e300;
                }
            }
            new_I_priors[number_of_regimes] = calculate_process_full_log_I_prior(i, m_particle.get_add_changepoint_index(), static_cast< unsigned int >(number_of_regimes)) + log(1 - m_particle.get_rho());
            full_log_acceptance_probability += calculate_log_of_sum_Bq_ratio(new_I_priors, old_I_prior);
        }
        // check if the old_I_prior_total equals the recorded value of full_I_prior
        //cout << "full I prior comparison " << old_I_prior_total - m_particle.get_log_full_I_prior() << endl;;
        //if (old_I_prior_total - m_particle.get_log_full_I_prior() > 0.0000001 || old_I_prior_total - m_particle.get_log_full_I_prior() < -0.0000001) {
        //    cout << "DD" << endl;
        //}
    }
    else if (u1 < m_d_k) {
        full_log_acceptance_probability += log(1 - m_p) - log(m_p);
        full_log_acceptance_probability += log(static_cast< double >(m_dimension)) - log(static_cast< double >(m_end - m_dimension + 1));
        full_log_acceptance_probability += m_dimension == 1 ? log(4.0) : 0;
        for (unsigned int i = 0; i < m_number_of_processes; i++) {
            double old_I_prior = calculate_process_full_log_I_prior_without_cp(i, m_h);
            size_t number_of_regimes = m_particle.get_number_of_regimes(i) - (m_particle.removing_full_changepoint_leaves_highest_regime_unobserved(i, m_particle.get_full_I_i_j(m_h, i)) ? 1 : 0);
            vector< double > new_I_priors = vector< double >(number_of_regimes + 1);
            for (unsigned int reg = 0; reg < number_of_regimes; reg++) {
                new_I_priors[reg] = calculate_process_full_log_I_prior_without_cp(i, m_h, reg);
                if ((reg == number_of_regimes - 1) && (m_particle.is_regime_unobserved(i, number_of_regimes - 1))) {
                    new_I_priors[reg] = -1e300;
                }
            }
            new_I_priors[number_of_regimes] = calculate_process_full_log_I_prior_without_cp(i, m_h, static_cast< unsigned int >(number_of_regimes)) + log(1 - m_particle.get_rho());
            full_log_acceptance_probability -= calculate_log_of_sum_Bq_ratio(new_I_priors, old_I_prior);
        }
    }
    else if (u1 < m_m_k) {
        for (unsigned int i = 0; i < m_number_of_processes; i++) {
            double old_I_prior = calculate_process_full_log_I_prior_without_cp(i, m_h);
            size_t number_of_regimes = m_particle.get_number_of_regimes(i) - (m_particle.removing_full_changepoint_leaves_highest_regime_unobserved(i, m_particle.get_full_I_i_j(m_h, i)) ? 1 : 0);
            vector< double > new_I_priors = vector< double >(number_of_regimes + 1);
            for (unsigned int reg = 0; reg < number_of_regimes; reg++) {
                new_I_priors[reg] = calculate_process_full_log_I_prior_without_cp(i, m_h, reg);
                if ((reg == number_of_regimes - 1) && (m_particle.is_regime_unobserved(i, number_of_regimes - 1))) {
                    new_I_priors[reg] = -1e300;
                }
            }
            new_I_priors[number_of_regimes] = calculate_process_full_log_I_prior_without_cp(i, m_h, static_cast< unsigned int >(number_of_regimes)) + log(1 - m_particle.get_rho());
            full_log_acceptance_probability += calculate_log_of_sum_Bq_ratio(new_I_priors, old_I_prior) - calculate_log_of_sum_Bq_ratio(new_I_priors, old_I_prior);
        }
    }
    else if (u1 < m_r_k) {
        full_log_acceptance_probability += 0;
    }
    else if (u1 < m_au_k) {
        full_log_acceptance_probability += 0;
    }
    return full_log_acceptance_probability;
}

double rj::calculate_process_full_log_I_prior(const double & process, const unsigned int & adding_index, const unsigned int & proposed_regime) {
    vector< unsigned int > regime_vector = vector< unsigned int >(0);
    for (unsigned int cp_idx = 0; cp_idx < m_dimension; cp_idx++) {
        regime_vector.push_back(m_particle.get_full_I_i_j(cp_idx, process));
    }
    size_t number_of_regimes = m_particle.get_number_of_regimes(process);
    // if we are calculating the full_log_I_prior for the process with a new regime added then add the new regime in
    if (adding_index != numeric_limits< unsigned int >::max()) {
        regime_vector.insert(regime_vector.begin() + adding_index, proposed_regime);
        if (proposed_regime == number_of_regimes) {
            number_of_regimes++;
        }
    }
    return calculate_full_log_I_prior(regime_vector, number_of_regimes);
}

double rj::calculate_process_full_log_I_prior_without_cp(const double & process, const unsigned int & removing_index, const unsigned int & proposed_regime) {
    vector< unsigned int > regime_vector = vector< unsigned int >(0);
    for (unsigned int cp_idx = 0; cp_idx < m_h; cp_idx++) {
        regime_vector.push_back(m_particle.get_full_I_i_j(cp_idx, process));
    }
    size_t number_of_regimes = m_particle.get_number_of_regimes(process);
    if (m_particle.removing_full_changepoint_leaves_highest_regime_unobserved(process, m_particle.get_full_I_i_j(m_h, process))) {
        number_of_regimes--;
    }
    if (proposed_regime != numeric_limits< unsigned int >::max()) {
        regime_vector.push_back(proposed_regime);
        if (proposed_regime == number_of_regimes) {
            number_of_regimes++;
        }
    }
    for (unsigned int cp_idx = m_h + 1; cp_idx < m_dimension; cp_idx++) {
        regime_vector.push_back(m_particle.get_full_I_i_j(cp_idx, process));
    }
    
    return calculate_full_log_I_prior(regime_vector, number_of_regimes);
}

double rj::calculate_full_log_I_prior(const vector< unsigned int > & regime_vector, const size_t & number_of_regimes) {
    vector< vector< double > > regime_transition_mat = generate_regime_transition_matrix(regime_vector, number_of_regimes);
    double num_reg = static_cast< double >(number_of_regimes), delta = m_particle.get_dirichlet_alpha();
    double full_log_I_prior = num_reg * (gsl_sf_lngamma(num_reg * delta) - num_reg * gsl_sf_lngamma(delta));
    for (unsigned int reg_0 = 0; reg_0 < number_of_regimes; reg_0++) {
        double iota_reg_0 = 0;
        for (unsigned int reg_1 = 0; reg_1 < number_of_regimes; reg_1++) {
            double iota_reg_0_reg_1 = regime_transition_mat[reg_0][reg_1];
            full_log_I_prior += gsl_sf_lngamma(delta + iota_reg_0_reg_1);
            iota_reg_0 += iota_reg_0_reg_1;
        }
        full_log_I_prior -= gsl_sf_lngamma(num_reg * delta + iota_reg_0);
    }
    return full_log_I_prior;
}

vector< vector< double > > rj::generate_regime_transition_matrix(const vector< unsigned int > & regime_vector, const size_t & number_of_regimes) {
    vector< vector< double > > regime_transition_matrix = vector< vector< double > >(number_of_regimes, vector< double >(number_of_regimes, 0.0));
    if (regime_vector.size() > 0) {
        regime_transition_matrix[0][regime_vector[0]]++;
        for (unsigned int idx = 0; idx < regime_vector.size() - 1; idx++) {
            regime_transition_matrix[regime_vector[idx]][regime_vector[idx + 1]]++;
        }
    }
    return regime_transition_matrix;
}

double rj::calculate_log_of_sum_Bq_ratio(const vector< double > & new_I_priors, const double & old_I_prior) {
    vector< double > a = new_I_priors;
    double log_of_sum_Bq = 0, max_a = a[0], temp_sum = 0;
    // subtract old_I_prior from each element and find the largest element
    for (unsigned int i = 0; i < a.size(); i++) {
        if (a[i] > max_a) {
            max_a = a[i];
        }
    }
    log_of_sum_Bq += max_a - old_I_prior;
    for (unsigned int i = 0; i < a.size(); i++) {
        temp_sum += exp(a[i] - max_a);
    }
    log_of_sum_Bq += log(temp_sum);
    return log_of_sum_Bq;
}

void rj::adding_basic_changepoint_setup(const unsigned int & trace_index) {
    unsigned long int lower_position_bound = ((trace_index == 0) ? 0 : m_separators[trace_index - 1]);
    unsigned long int upper_position_bound = ((trace_index == m_number_of_traces - 1) ? m_end + 1 : m_separators[trace_index]);
    m_new_changepoint_position = gsl_rng_uniform_int(r, upper_position_bound - lower_position_bound - 1) + lower_position_bound + 1;
    while (m_particle.does_changepoint_exist_in_particle(m_new_changepoint_position)) {
        m_new_changepoint_position = gsl_rng_uniform_int(r, upper_position_bound - lower_position_bound - 1) + lower_position_bound + 1;
    }
    m_log_proposal_ratio = m_particle.calculate_and_get_add_cp_proposal_ratio(upper_position_bound - lower_position_bound - 1, trace_index, true);
    m_adding_changepoint = changepoint(m_new_changepoint_position);
    unsigned int add_cp_index = m_particle.get_add_changepoint_index();
    m_left_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(add_cp_index - 1).get_position(), m_new_changepoint_position);
    if (add_cp_index == m_dimension) {
        m_right_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_new_changepoint_position, m_end + 1);
    } else {
        m_right_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_new_changepoint_position, m_particle.get_changepoint(add_cp_index).get_position());
    }
    m_log_likelihood_ratio = m_left_log_likelihood + m_right_log_likelihood - m_particle.get_changepoint(add_cp_index - 1).get_log_likelihood();
    m_log_k_prior_ratio = m_particle.calculate_and_get_add_cp_k_prior_ratio();
    m_log_acceptance_prob = m_log_proposal_ratio + m_log_k_prior_ratio + m_log_likelihood_ratio;
}

void rj::removing_basic_changepoint_setup(const unsigned int & trace_index) {
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    m_log_proposal_ratio = m_particle.calculate_and_get_remove_cp_proposal_ratio(trace_dimension, trace_index, true);
    if (m_h == m_dimension - 1) {
        m_merged_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(m_h - 1).get_position(), m_end + 1);
    } else {
        m_merged_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_particle.get_changepoint(m_h - 1).get_position(), m_particle.get_changepoint(m_h + 1).get_position());
    }
    m_log_likelihood_ratio = m_merged_log_likelihood - m_particle.get_changepoint(m_h - 1).get_log_likelihood() - m_particle.get_changepoint(m_h).get_log_likelihood();
    m_log_k_prior_ratio = m_particle.calculate_and_get_remove_cp_k_prior_ratio();
    m_log_acceptance_prob = m_log_proposal_ratio + m_log_k_prior_ratio + m_log_likelihood_ratio;
}

void rj::moving_basic_changepoint_setup(const unsigned int & trace_index) {
    // choose which changepoint to move (we know that there is at least one that can be moved from our definitions of b_k, d_k, etc.)
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    unsigned long int tau_h_minus_1 = m_particle.get_changepoint(m_h - 1).get_position();
    unsigned long int tau_h_plus_1;
    if (m_h == m_dimension - 1) {
        tau_h_plus_1 = m_end + 1;
    } else {
        tau_h_plus_1 = m_particle.get_changepoint(m_h + 1).get_position();
    }
    m_new_changepoint_position = gsl_rng_uniform_int(r, tau_h_plus_1 - tau_h_minus_1 - 1) + tau_h_minus_1 + 1;
    if (tau_h_plus_1 - tau_h_minus_1 > 2) { //if there is somewhere to move the changepoint to. E.g. 16, 17, 18 move 17 - there is nowhere to move it to!
        while(m_particle.does_changepoint_exist_in_particle(m_new_changepoint_position, m_h, m_h)) {
            m_new_changepoint_position = gsl_rng_uniform_int(r, tau_h_plus_1 - tau_h_minus_1 - 1) + tau_h_minus_1 + 1;
        }
    }
    m_log_proposal_ratio = 0;
    m_log_k_prior_ratio = 0;
    m_left_log_likelihood = m_pm_ptr->calculate_log_likelihood(tau_h_minus_1, m_new_changepoint_position);
    m_right_log_likelihood = m_pm_ptr->calculate_log_likelihood(m_new_changepoint_position, tau_h_plus_1);
    m_log_likelihood_ratio = m_left_log_likelihood + m_right_log_likelihood - m_particle.get_changepoint(m_h - 1).get_log_likelihood() - m_particle.get_changepoint(m_h).get_log_likelihood();
    m_log_acceptance_prob = m_log_likelihood_ratio;
}

void rj::basic_acceptance_procedure(const double & u1) {
	m_particle.increase_log_likelihood(m_log_likelihood_ratio);
	m_particle.increase_log_k_prior(m_log_k_prior_ratio);
	if (u1 < m_b_k) {
		m_adding_changepoint.set_log_likelihood(m_right_log_likelihood);
		m_particle.add_changepoint(m_particle.get_add_changepoint_index(), m_adding_changepoint);
		m_particle.get_changepoint(m_particle.get_add_changepoint_index() - 1).set_log_likelihood(m_left_log_likelihood);
	}
	else if (u1 < m_d_k) {
		m_particle.remove_changepoint(m_h);
		m_particle.get_changepoint(m_h - 1).set_log_likelihood(m_merged_log_likelihood);
	}
	else if (u1 < m_m_k) {
        m_particle.move_changepoint(m_h, m_new_changepoint_position, m_left_log_likelihood, m_right_log_likelihood);
	}
}

void rj::basic_recording_procedure() {
	m_recorded_basic_dimensions.push_back(m_dimension);
	vector< unsigned long int > changepoint_hist = m_particle.calculate_and_get_changepoint_histogram(m_number_of_changepoint_bins);
	size_t size_of_recorded_changepoints = m_recorded_basic_changepoints.size();
	if (size_of_recorded_changepoints > 0) { //have we started recording changepoint histograms, or is this the first time?
		for (unsigned long int i = 0; i < size_of_recorded_changepoints; i++) {
			m_recorded_basic_changepoints[i] += changepoint_hist[i];
		}
	} else {
		m_recorded_basic_changepoints = changepoint_hist;
	}
    m_recorded_basic_log_posteriors.push_back(m_particle.get_basic_log_posterior());
}

void rj::run_basic_simulation(){
    m_particle.print_likelihood();
    cout << "the posterior is " << m_particle.get_basic_log_posterior() << endl << endl;
    //if (abs(m_particle.calculate_and_get_basic_log_posterior(m_number_of_processes) - m_particle.get_basic_log_posterior()) > 0.0001) {
    //    cout << m_particle.calculate_and_get_basic_log_posterior(m_number_of_processes) - m_particle.get_basic_log_posterior() << endl;
    //}
    cout << "dimension: " << m_particle.get_dimension() << endl;
    double accepted_adds = 0, accepted_deaths = 0, accepted_moves = 0;
    double attempted_adds = 0, attempted_deaths = 0, attempted_moves = 0;
    for (unsigned long int iteration = 0; iteration < m_basic_burnin; iteration++) {
        m_dimension = m_particle.get_dimension();
        unsigned int trace_index = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_number_of_traces));
        unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
        if (trace_dimension > 0) {
            m_b_k = 1.0 / 3.0, m_d_k = 2.0 / 3.0, m_m_k = 1.0;
        } else {
            m_b_k = 1.0;
        }
        
        /*if (iteration == 100 || iteration == 250 || iteration == 500 || iteration == 750 || iteration == 1000 || iteration == 10000) {
            cout << "iteration: " << iteration << endl;
            cout << "dimension: " << m_dimension << endl;
            cout << "birth: " << accepted_adds << " / " << attempted_adds << " = " << accepted_adds / attempted_adds << endl;
            cout << "death: " << accepted_deaths << " / " << attempted_deaths << " = " << accepted_deaths / attempted_deaths << endl;
            cout << "move: " << accepted_moves << " / " << attempted_moves << " = " << accepted_moves / attempted_moves << endl;
        }*/
        
        //m_particle.check_separator_changepoints();

        double u1 = gsl_ran_flat(r, 0, 1);
        if (u1 < m_b_k) { //birth
            adding_basic_changepoint_setup(trace_index);
            attempted_adds += 1;
        } else if (u1 < m_d_k) { //death
            removing_basic_changepoint_setup(trace_index);
            attempted_deaths += 1;
        } else if (u1 < m_m_k) { //move
            moving_basic_changepoint_setup(trace_index);
            attempted_moves += 1;
        }
        
        if (m_log_acceptance_prob > 0 || (log(gsl_ran_flat(r, 0, 1)) < m_log_acceptance_prob)) {
			basic_acceptance_procedure(u1);
            if (u1 < m_b_k) { //birth
                accepted_adds += 1;
            } else if (u1 < m_d_k) { //death
                accepted_deaths += 1;
            } else if (u1 < m_m_k) { //move
                accepted_moves += 1;
            }
        }
        
        //if (abs(m_particle.calculate_and_get_basic_log_posterior(m_number_of_processes) - m_particle.get_basic_log_posterior()) > 0.0001) {
        //    cout << iteration << '\t' << m_particle.calculate_and_get_basic_log_posterior(m_number_of_processes) - m_particle.get_basic_log_posterior() << endl;
        //}
    }
    
    for (unsigned long int iteration = 0; iteration < m_basic_iterations; iteration++) {
        m_dimension = m_particle.get_dimension();
        unsigned int trace_index = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_number_of_traces));
        unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
        if (trace_dimension > 0) {
            m_b_k = 1.0 / 3.0, m_d_k = 2.0 / 3.0, m_m_k = 1.0;
        } else {
            m_b_k = 1.0;
        }
		
		double u1 = gsl_ran_flat(r, 0, 1);
		if (u1 < m_b_k) { //birth
			adding_basic_changepoint_setup(trace_index);
		}
		else if (u1 < m_d_k) { //death
			removing_basic_changepoint_setup(trace_index);
		}
		else if (u1 < m_m_k) { //move
			moving_basic_changepoint_setup(trace_index);
		}

		if (m_log_acceptance_prob > 0 || (log(gsl_ran_flat(r, 0, 1)) < m_log_acceptance_prob)) {
			basic_acceptance_procedure(u1);
		}
        
        if (m_recording_basic_samples) {
            double log_posterior = m_particle.get_basic_log_posterior();
            if (m_basic_MAP_log_posterior < log_posterior) {
                m_basic_MAP_particle = m_particle;
                m_basic_MAP_dimension = m_particle.get_dimension();
                m_basic_MAP_log_posterior = log_posterior;
            }
        }

		if (m_recording_basic_samples && (iteration % m_basic_thinning == 0)) { //store sample
			basic_recording_procedure();
		}
    }
    m_dimension = m_particle.get_dimension();
    cout << "ending basic changepoint stage" << endl;
    m_particle.print_likelihood();
    cout << "the posterior is " << m_particle.get_basic_log_posterior() << endl << endl;
}

void rj::write_basic_MAP_changepoints_to_file(const string & MAP_cps_Filename){
    ofstream OutputStream(MAP_cps_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (unsigned int cp_index = 0; cp_index < m_basic_MAP_dimension; cp_index++) {
        OutputStream << m_basic_MAP_particle.get_changepoint(cp_index).get_position() << "\n";
    }
    OutputStream.close();
}

void rj::write_binary_MAP_changepoints_to_file(const string & MAP_cps_Filename){
    ofstream OutputStream(MAP_cps_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (unsigned int cp_index = 0; cp_index < m_binary_MAP_dimension; cp_index++) {
        OutputStream << m_binary_MAP_particle.get_changepoint(cp_index).get_position() * m_pm_ptr->get_diff() << "\t";
        for (unsigned int process = 0; process < m_number_of_processes; process++) {
            OutputStream << m_binary_MAP_particle.get_binary_I_i_j(cp_index, process) << "\t";
        }
        OutputStream << endl;
    }
    OutputStream.close();
}

void rj::write_full_MAP_changepoints_to_file(const string & MAP_cps_Filename) {
    ofstream OutputStream(MAP_cps_Filename, ios::out);
    OutputStream.precision(10);
    OutputStream << m_full_MAP_log_posterior << endl;
    OutputStream << m_full_MAP_dimension << endl;
    for (unsigned int cp_index = 0; cp_index < m_full_MAP_dimension; cp_index++) {
        OutputStream << m_full_MAP_particle.get_changepoint(cp_index).get_position() * m_pm_ptr->get_diff() << "\t";
        for (unsigned int process = 0; process < m_number_of_processes; process++) {
            OutputStream << m_full_MAP_particle.get_full_I_i_j(cp_index, process) << "\t";
        }
        OutputStream << endl;
    }
    OutputStream.close();
}

void rj::write_basic_dimension_distribution_to_file( const string & dimension_distribution_Filename ){
    ofstream OutputStream( dimension_distribution_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    unsigned long int max_dim = 0;
    for (size_t i = 0; i < m_recorded_basic_dimensions.size(); i++) {
        if (max_dim < m_recorded_basic_dimensions[i]){
            max_dim = m_recorded_basic_dimensions[i];
        }
    }
    vector< double > dimension_counts(max_dim + 1, 0.0);
    for (size_t i = 0; i < m_recorded_basic_dimensions.size(); i++) {
        dimension_counts[m_recorded_basic_dimensions[i]]++;
    }
    for (size_t j = 0; j <= max_dim; j++) {
        dimension_counts[j] /= static_cast< double >(m_recorded_basic_dimensions.size());
    }
    for (unsigned int k = 0; k <= max_dim; k++) {
        OutputStream << k << "\t" << dimension_counts[k] << "\n";
    }
    OutputStream.close();
}

void rj::write_binary_dimension_distribution_to_file(const string & dimension_distribution_Filename) {
    ofstream OutputStream( dimension_distribution_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    // calculate the maximum observed dimension
    unsigned long int max_dim = 0;
    for (size_t i = 0; i < m_recorded_binary_dimensions.size(); i++) {
        if (max_dim < m_recorded_binary_dimensions[i]) {
            max_dim = m_recorded_binary_dimensions[i];
        }
    }
	// calculate a histogram of dimension counts
    vector< double > dimension_counts(max_dim + 1, 0.0);
    for (size_t i = 0; i < m_recorded_binary_dimensions.size(); i++) {
        dimension_counts[m_recorded_binary_dimensions[i]]++;
    }
    // convert from a frequency to a probability
    for (size_t j = 0; j <= max_dim; j++) {
        dimension_counts[j] /= static_cast< double >(m_recorded_binary_dimensions.size());
    }
    for (unsigned int k = 0; k <= max_dim; k++) {
        OutputStream << k << "\t" << dimension_counts[k] << "\n";
    }
    OutputStream.close();
}

void rj::write_full_dimension_distribution_to_file(const string & dimension_distribution_Filename){
    ofstream OutputStream(dimension_distribution_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    // calculate what the largest observed dimension was
    unsigned long int max_dim = 0;
    for (size_t i = 0; i < m_recorded_full_dimensions.size(); i++){
        if (max_dim < m_recorded_full_dimensions[i]){
            max_dim = m_recorded_full_dimensions[i];
        }
    }
	// create a histogram of dimension counts
    vector< double > dimension_counts(max_dim + 1, 0.0);
    for (unsigned int i = 0; i < m_recorded_full_dimensions.size(); i++){
        dimension_counts[m_recorded_full_dimensions[i]]++;
    }
    // convert the counts from a frequency to a probability
    for (size_t j = 0; j <= max_dim; j++){
        dimension_counts[j] /= static_cast< double >(m_recorded_full_dimensions.size());
    }
    for( unsigned int k = 0; k <= max_dim; k++ ){
        OutputStream << k << "\t" << dimension_counts[k] << "\n";
    }
    OutputStream.close();
}

void rj::write_full_effective_dimension_distribution_to_file(const string & effective_dimension_distribution_Filename){
    ofstream OutputStream(effective_dimension_distribution_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    // calculate what the largest observed dimension was
    unsigned long int max_dim = 0;
    for (size_t i = 0; i < m_recorded_full_effective_dimensions.size(); i++){
        if (max_dim < m_recorded_full_effective_dimensions[i]){
            max_dim = m_recorded_full_effective_dimensions[i];
        }
    }
    // create a histogram of dimension counts
    vector< double > dimension_counts(max_dim + 1, 0.0);
    for (unsigned int i = 0; i < m_recorded_full_effective_dimensions.size(); i++){
        dimension_counts[m_recorded_full_effective_dimensions[i]]++;
    }
    // convert the counts from a frequency to a probability
    for (size_t j = 0; j <= max_dim; j++){
        dimension_counts[j] /= static_cast< double >(m_recorded_full_effective_dimensions.size());
    }
    for( unsigned int k = 0; k <= max_dim; k++ ){
        OutputStream << k << "\t" << dimension_counts[k] << "\n";
    }
    OutputStream.close();
}

void rj::write_basic_changepoints_distribution_to_file(const string & changepoints_distribution_Filename){
    ofstream OutputStream( changepoints_distribution_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (size_t i = 0; i < m_recorded_basic_changepoints.size(); i++){
        OutputStream << static_cast< double >(m_recorded_basic_changepoints[i]) / static_cast< double >(m_basic_iterations / m_basic_thinning) << "\n";
    }
}

void rj::write_binary_changepoints_distribution_to_file(const string & changepoints_distribution_Filename){
    ofstream OutputStream( changepoints_distribution_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    if (m_recorded_binary_changepoints.size() != m_number_of_processes * m_number_of_changepoint_bins){
        cerr << "don't match";
    }
    for (unsigned int j = 0; j < m_number_of_processes; j++ ){
        for (unsigned int index = 0; index < m_number_of_changepoint_bins; index++) {
            OutputStream << static_cast< double >(m_recorded_binary_changepoints[m_number_of_changepoint_bins * j + index]) / static_cast< double >(m_binary_iterations / m_binary_thinning) << "\t";
        }
        OutputStream << endl;
    }
}

void rj::write_full_changepoints_distribution_to_file(const string & changepoints_distribution_Filename, const unsigned long int number_of_iterations){
    ofstream OutputStream(changepoints_distribution_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(5);
    if (m_recorded_full_changepoints.size() != m_number_of_processes * m_number_of_changepoint_bins){
        cerr << "don't match";
    }
    for (unsigned int j = 0; j < m_number_of_processes; j++ ){
        for (unsigned int index = 0; index < m_number_of_changepoint_bins; index++) {
            OutputStream << static_cast< double >(m_recorded_full_changepoints[m_number_of_changepoint_bins * j + index]) / static_cast< double >(number_of_iterations / m_full_thinning) << "\t";
        }
        OutputStream << endl;
    }
}

void rj::write_number_of_regimes_to_file(const string & number_of_regimes_Filename) {
    ofstream OutputStream(number_of_regimes_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (unsigned int proc = 0; proc < m_number_of_processes; proc++) {
        // calculate the largest recorded number of regimes for this process
        unsigned long int max_r = 1;
        for (size_t i = 0; i < m_recorded_number_of_regimes.size(); i++){
            if (max_r < m_recorded_number_of_regimes[i][proc]){
                max_r = m_recorded_number_of_regimes[i][proc];
            }
        }
        // create a histogram of regime counts
        vector< double > regime_counts(max_r + 1, 0.0);
        for (unsigned int i = 0; i < m_recorded_number_of_regimes.size(); i++){
            regime_counts[static_cast< double >(m_recorded_number_of_regimes[i][proc])]++;
        }
        // convert the counts from a frequency to a probability
        for (size_t j = 0; j <= max_r; j++){
            regime_counts[j] /= static_cast< double >(m_recorded_full_dimensions.size());
        }
        for( unsigned int k = 1; k <= max_r; k++ ){
            OutputStream << proc << "\t" << k << "\t" << regime_counts[k] << "\n";
        }
    }
    OutputStream.close();
}

void rj::write_number_of_observed_regimes_to_file(const string & number_of_observed_regimes_Filename) {
    ofstream OutputStream(number_of_observed_regimes_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (unsigned int proc = 0; proc < m_number_of_processes; proc++) {
        // calculate the largest recorded number of regimes for this process
        unsigned long int max_r = 1;
        for (size_t i = 0; i < m_recorded_number_of_observed_regimes.size(); i++){
            if (max_r < m_recorded_number_of_observed_regimes[i][proc]){
                max_r = m_recorded_number_of_observed_regimes[i][proc];
            }
        }
        // create a histogram of regime counts
        vector< double > regime_counts(max_r + 1, 0.0);
        for (unsigned int i = 0; i < m_recorded_number_of_observed_regimes.size(); i++){
            regime_counts[static_cast< double >(m_recorded_number_of_observed_regimes[i][proc])]++;
        }
        // convert the counts from a frequency to a probability
        for (size_t j = 0; j <= max_r; j++){
            regime_counts[j] /= static_cast< double >(m_recorded_full_dimensions.size());
        }
        for( unsigned int k = 0; k <= max_r; k++ ){
            OutputStream << proc << "\t" << k << "\t" << regime_counts[k] << "\n";
        }
    }
    OutputStream.close();
}

void rj::write_basic_log_posterior_trace_to_file(const string & log_posterior_trace_Filename){
    ofstream OutputStream(log_posterior_trace_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (size_t i = 0; i < m_recorded_basic_log_posteriors.size(); i++ ){
        OutputStream << m_recorded_basic_log_posteriors[i] << "\n";
    }
}

void rj::write_binary_log_posterior_trace_to_file(const string & log_posterior_trace_Filename){
    ofstream OutputStream(log_posterior_trace_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (size_t i = 0; i < m_recorded_binary_log_posteriors.size(); i++ ){
        OutputStream << m_recorded_binary_log_posteriors[i] << "\n";
    }
}

void rj::write_full_log_posterior_trace_to_file(const string & log_posterior_trace_Filename){
    ofstream OutputStream(log_posterior_trace_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (size_t i = 0; i < m_recorded_full_log_posteriors.size(); i++ ){
        OutputStream << m_recorded_full_log_posteriors[i] << "\n";
    }
}

void rj::write_basic_dimension_trace_to_file(const string & dimension_trace_Filename) {
    ofstream OutputStream(dimension_trace_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (size_t i = 0; i < m_recorded_basic_dimensions.size(); i++ ){
        OutputStream << m_recorded_basic_dimensions[i] << "\n";
    }
}

void rj::write_dimension_trace_to_file(const string & dimension_trace_Filename) {
    ofstream OutputStream(dimension_trace_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (size_t i = 0; i < m_recorded_full_dimensions.size(); i++ ){
        OutputStream << m_recorded_full_dimensions[i] << "\n";
    }
}

void rj::write_number_of_regimes_trace_to_file(const string & number_of_regimes_trace_Filename) {
    ofstream OutputStream(number_of_regimes_trace_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    for (size_t i = 0; i < m_recorded_number_of_regimes.size(); i++ ) {
        for (size_t j = 0; j < m_recorded_number_of_regimes[i].size() - 1; j++) {
            OutputStream << m_recorded_number_of_regimes[i][j] << ",";
        }
        OutputStream << m_recorded_number_of_regimes[i][m_recorded_number_of_regimes[i].size() - 1] << endl;
    }
}

void rj::write_full_acceptance_probabilities_to_file(const string & acceptance_probabilities_Filename) {
    ofstream OutputStream(acceptance_probabilities_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(10);
    OutputStream << "full birth acceptance probability" << "\t" << m_recorded_full_birth_acceptances / m_recorded_full_birth_proposals << " = " << static_cast< unsigned int >(m_recorded_full_birth_acceptances) << " / " << static_cast< unsigned int >(m_recorded_full_birth_proposals) << endl;
    OutputStream << "full death acceptance probability" << "\t" << m_recorded_full_death_acceptances / m_recorded_full_death_proposals << " = " << static_cast< unsigned int >(m_recorded_full_death_acceptances) << " / " << static_cast< unsigned int >(m_recorded_full_death_proposals) << endl;
    OutputStream << "full move acceptance probability" << "\t" << m_recorded_full_move_acceptances / m_recorded_full_move_proposals << " = " << static_cast< unsigned int >(m_recorded_full_move_acceptances) << " / " << static_cast< unsigned int >(m_recorded_full_move_proposals) << endl;
    OutputStream << "full resample acceptance probability" << "\t" << m_recorded_full_resample_acceptances / m_recorded_full_resample_proposals << " = " << static_cast< unsigned int >(m_recorded_full_resample_acceptances) << " / " << static_cast< unsigned int >(m_recorded_full_resample_proposals) << endl;
    OutputStream << "full alter unobserveds acceptance probability" << "\t" << m_recorded_full_unobserveds_acceptances / m_recorded_full_unobserveds_proposals << " = " << static_cast< unsigned int >(m_recorded_full_unobserveds_acceptances) << " / " << static_cast< unsigned int >(m_recorded_full_unobserveds_proposals) << endl;
    OutputStream << "full overall acceptance probability" << "\t" << (m_recorded_full_birth_acceptances + m_recorded_full_death_acceptances + m_recorded_full_move_acceptances + m_recorded_full_resample_acceptances + m_recorded_full_unobserveds_acceptances) / (m_recorded_full_birth_proposals + m_recorded_full_death_proposals + m_recorded_full_move_proposals + m_recorded_full_resample_proposals + m_recorded_full_unobserveds_proposals) << " = " << static_cast< unsigned int >(m_recorded_full_birth_acceptances + m_recorded_full_death_acceptances + m_recorded_full_move_acceptances + m_recorded_full_resample_acceptances + m_recorded_full_unobserveds_acceptances) << " / " << static_cast< unsigned int >(m_recorded_full_birth_proposals + m_recorded_full_death_proposals + m_recorded_full_move_proposals + m_recorded_full_resample_proposals + m_recorded_full_unobserveds_proposals);
}

void rj::write_similarity_matrix_to_file(const string & similarity_matrix_Filename) {
    ofstream OutputStream(similarity_matrix_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(5);
    if (m_recorded_similarity_matrix.size() != m_number_of_traces) {
        cerr << "similarity matrix size doesn't match";
    }
    for (unsigned int trace_0 = 0; trace_0 < m_number_of_traces; trace_0++) {
        for (unsigned int trace_1 = 0; trace_1 < m_number_of_traces; trace_1++) {
            if (trace_0 == trace_1) {
                OutputStream << 1 << '\t';
            }
            else {
                OutputStream << m_recorded_similarity_matrix[trace_0][trace_1] / (static_cast< double >(m_full_iterations / m_full_thinning) * static_cast< double >(m_observations_in_each_trace[trace_0] + m_observations_in_each_trace[trace_1])) << '\t';
            }
        }
        OutputStream << endl;
    }
}

void rj::write_min_proportion_similarity_matrix_to_file(const string & min_proportion_similarity_matrix_Filename) {
    ofstream OutputStream(min_proportion_similarity_matrix_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(5);
    if (m_recorded_min_proportion_similarity_matrix.size() != m_number_of_traces) {
        cerr << "similarity matrix size doesn't match";
    }
    for (unsigned int trace_0 = 0; trace_0 < m_number_of_traces; trace_0++) {
        for (unsigned int trace_1 = 0; trace_1 < m_number_of_traces; trace_1++) {
            if (trace_0 == trace_1) {
                OutputStream << 1 << '\t';
            }
            else {
                OutputStream << m_recorded_min_proportion_similarity_matrix[trace_0][trace_1] / static_cast< double >(m_full_iterations / m_full_thinning) << '\t';
            }
        }
        OutputStream << endl;
    }
}

void rj::write_similarity_matrices_to_file(const string & similarity_matrices_Filename) {
    ofstream OutputStream(similarity_matrices_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(5);
	// start by making a header for the file, detailing which pair of traces we are recording the similarity of
	for (unsigned int trace_0 = 0; trace_0 < m_number_of_traces; trace_0++) {
		for (unsigned int trace_1 = trace_0 + 1; trace_1 < m_number_of_traces; trace_1++) {
			if (trace_0 == m_number_of_traces - 2 && trace_1 == m_number_of_traces - 1) {
				OutputStream << "traces_" << trace_0 << "_" << trace_1 << endl;
			}
			else {
				OutputStream << "traces_" << trace_0 << "_" << trace_1 << '\t';
			}
		}
	}

	for (unsigned iter = 0; iter < m_recorded_similarity_matrices.size(); iter++) {
		for (unsigned int trace_0 = 0; trace_0 < m_number_of_traces; trace_0++) {
			for (unsigned int trace_1 = trace_0 + 1; trace_1 < m_number_of_traces; trace_1++) {
				if (trace_0 == m_number_of_traces - 2 && trace_1 == m_number_of_traces - 1) {
					OutputStream << m_recorded_similarity_matrices[iter][trace_0][trace_1] / static_cast< double >(m_observations_in_each_trace[trace_0] + m_observations_in_each_trace[trace_1]) << endl;
				}
				else {
					OutputStream << m_recorded_similarity_matrices[iter][trace_0][trace_1] / static_cast< double >(m_observations_in_each_trace[trace_0] + m_observations_in_each_trace[trace_1]) << '\t';
				}
			}
		}
	}
}

void rj::write_min_proportion_similarity_matrices_to_file(const string & min_proportion_similarity_matrices_Filename) {
	ofstream OutputStream(min_proportion_similarity_matrices_Filename, ios::out);
	OutputStream << setiosflags(ios::fixed);
	OutputStream.precision(5);
	// start by making a header for the file, detailing which pair of traces we are recording the similarity of
	for (unsigned int trace_0 = 0; trace_0 < m_number_of_traces; trace_0++) {
		for (unsigned int trace_1 = trace_0 + 1; trace_1 < m_number_of_traces; trace_1++) {
			if (trace_0 == m_number_of_traces - 2 && trace_1 == m_number_of_traces - 1) {
				OutputStream << "traces_" << trace_0 << "_" << trace_1 << endl;
			}
			else {
				OutputStream << "traces_" << trace_0 << "_" << trace_1 << '\t';
			}
		}
	}

	for (unsigned iter = 0; iter < m_recorded_similarity_matrices.size(); iter++) {
		for (unsigned int trace_0 = 0; trace_0 < m_number_of_traces; trace_0++) {
			for (unsigned int trace_1 = trace_0 + 1; trace_1 < m_number_of_traces; trace_1++) {
				if (trace_0 == m_number_of_traces - 2 && trace_1 == m_number_of_traces - 1) {
					OutputStream << m_recorded_min_proportion_similarity_matrices[iter][trace_0][trace_1] << endl;
				}
				else {
					OutputStream << m_recorded_min_proportion_similarity_matrices[iter][trace_0][trace_1] << '\t';
				}
			}
		}
	}
}

void rj::write_association_matrix_to_file(const string & association_matrix_Filename) {
    ofstream OutputStream(association_matrix_Filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    OutputStream.precision(5);
    
    // divide the values in the matrices by the number of records made
    double num_recordings = static_cast< double >(m_full_iterations / m_full_thinning);
    
    OutputStream << "x\ty\tSimilarity\tProcess";
    for (unsigned int process = 0; process < m_number_of_processes; process++) {
        for (unsigned long int ind_0 = 0; ind_0 < m_number_of_association_matrix_bins; ind_0++) {
            OutputStream << endl << ind_0 * m_end / m_number_of_association_matrix_bins << "\t" << ind_0 * m_end / m_number_of_association_matrix_bins << "\t" << m_association_matrices[process][ind_0][ind_0] / num_recordings << "\t" << process;
            for (unsigned long int ind_1 = ind_0 + 1; ind_1 < m_number_of_association_matrix_bins; ind_1++) {
                OutputStream << endl << ind_0 * m_end / m_number_of_association_matrix_bins << "\t" << ind_1 * m_end / m_number_of_association_matrix_bins << "\t" << m_association_matrices[process][ind_0][ind_1] / num_recordings << "\t" << process;
                OutputStream << endl << ind_1 * m_end / m_number_of_association_matrix_bins << "\t" << ind_0 * m_end / m_number_of_association_matrix_bins << "\t" << m_association_matrices[process][ind_0][ind_1] / num_recordings << "\t" << process;
            }
        }
    }
}

// work out the likelihoods for each segment for each process. Create the binaries and set the value of alpha for binary marked vector sampling
void rj::convert_basic_particle_to_binary_particle(const double & beta_alpha) {
    /*m_number_of_processes = m_pm_ptr->get_number_of_processes();
    vector< vector< double > > log_likelihoods;
    for (unsigned int proc = 0; proc < m_number_of_processes; proc++) {
        log_likelihoods.push_back(vector< double >(m_dimension + 1));
        if (m_dimension > 0) {
            log_likelihoods[proc][0] = m_pm_ptr->calculate_log_likelihood(proc, m_particle.get_changepoint(-1).get_position(), m_particle.get_changepoint(0).get_position());
            for (unsigned int i = 0; i < m_dimension - 1; i++) {
                log_likelihoods[proc][i + 1] = m_pm_ptr->calculate_log_likelihood(proc, m_particle.get_changepoint(i).get_position(), m_particle.get_changepoint(i + 1).get_position());
            }
            log_likelihoods[proc][m_dimension] = m_pm_ptr->calculate_log_likelihood(proc, m_particle.get_changepoint(m_dimension - 1).get_position(), m_end + 1);
        }
        else {
            log_likelihoods[proc][0] = m_pm_ptr->calculate_log_likelihood(proc, m_particle.get_changepoint(-1).get_position(), m_end + 1);
        }
    }*/
    m_number_of_processes = m_pm_ptr->get_number_of_processes();
	m_particle.initiate_binaries(m_number_of_processes);
	m_particle.set_log_binary_I_prior(0);
    m_particle.set_log_likelihood(0);
    m_particle.set_beta_alpha(beta_alpha);
	set_binary_marked_vectors();//m_number_of_processes, log_likelihoods);
    m_particle.calculate_and_set_log_binary_I_prior(m_number_of_processes);
    //cout << calculate_total_binary_log_likelihood() << endl;
    m_binary_left_log_likelihood = vector< double >(m_number_of_processes);
    m_binary_right_log_likelihood = vector< double >(m_number_of_processes);
}

void rj::set_binary_marked_vectors() {
	m_particle.set_all_binary_marked_vectors_equal_to_0_vectors(m_number_of_processes);
	for (unsigned int j = 0; j < m_number_of_processes; j++) {
		unsigned int trace_index = 0;
		// start with cp_index = -1
		int cp_index = -1;
		unsigned long int left_cp_position = m_particle.get_changepoint(cp_index).get_position();
        unsigned long int right_cp_position;
        if (m_particle.get_dimension() == 0) {
            right_cp_position = m_end + 1;
        }
        else {
            right_cp_position = m_particle.get_changepoint(cp_index + 1).get_position();
        }
		// we get the sufficient statistics and number of observations for all the processes here when we are only interested in the statistics for process j, but don't want to prematurely optimise here
        double log_likelihood = m_pm_ptr->calculate_log_likelihood(j, left_cp_position, right_cp_position);
        /*vector< vector< double > > stats_1(m_number_of_processes);
		vector< vector< double > > stats_2(m_number_of_processes);
		m_pm_ptr->get_cumulative_sufficient_data(left_cp_position, stats_1);
		m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, stats_2);
		for (size_t i = 0; i < stats_1.size(); i++) {
			stats_2[j][i] -= stats_1[j][i];
		}
		double number_of_observations;
		m_pm_ptr->get_number_of_observations(stats_2[j], number_of_observations);
		double log_likelihood = m_pm_ptr->calculate_log_likelihood(j, stats_2[j]);*/
		double log_prior = m_particle.calculate_and_get_adding_binary_log_I_prior_ratio(j, true, trace_index, cp_index)[1];

		// add this binary to the particle
        m_particle.increase_log_likelihood(log_likelihood);
        m_particle.increase_log_binary_I_prior(log_prior);
		m_particle.add_new_binary(j, -1, log_likelihood);
		m_particle.increase_log_binary_I_prior(log_prior);

		for (int cp_index = 0; cp_index < m_particle.get_dimension(); cp_index++) {
			unsigned long int left_cp_position = m_particle.get_changepoint(cp_index).get_position();
            unsigned long int right_cp_position;
            if (cp_index == m_particle.get_dimension() - 1) {
                right_cp_position = m_end + 1;
            }
            else {
                right_cp_position = m_particle.get_changepoint(cp_index + 1).get_position();
            }
			bool new_trace = false;
			if (m_particle.is_changepoint_index_separator_index(cp_index)) {
				trace_index++;
				new_trace = true;
			}
			
			// we get the sufficient statistics and number of observations for all the processes here when we are only interested in the statistics for process j, but don't want to prematurely optimise here
            double right_log_likelihood = m_pm_ptr->calculate_log_likelihood(j, left_cp_position, right_cp_position);
            /*vector< vector< double > > stats_1(m_number_of_processes);
			vector< vector< double > > stats_2(m_number_of_processes);
			m_pm_ptr->get_cumulative_sufficient_data(left_cp_position, stats_1);
			m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, stats_2);
			for (size_t i = 0; i < stats_1.size(); i++) {
				stats_2[j][i] -= stats_1[j][i];
			}
			double right_log_likelihood = m_pm_ptr->calculate_log_likelihood(j, stats_2[j]);*/

            unsigned long int previous_binary_left_position = m_particle.get_changepoint(m_particle.get_binary_left_index(j, cp_index)).get_position();
            double combined_log_likelihood = m_pm_ptr->calculate_log_likelihood(j, previous_binary_left_position, right_cp_position);
			/*vector< vector< double > > stats_3(m_number_of_processes);
			vector< vector< double > > stats_4(m_number_of_processes);
			m_pm_ptr->get_cumulative_sufficient_data(previous_binary_left_position, stats_3);
			m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, stats_4);
			for (size_t i = 0; i < stats_1.size(); i++) {
				stats_4[j][i] -= stats_3[j][i];
			}
			double combined_log_likelihood = m_pm_ptr->calculate_log_likelihood(j, stats_4[j]);*/

			vector< double > log_B(2);
			vector< double > log_q(2);
			log_B[0] = combined_log_likelihood - m_particle.get_binary_log_likelihood(j, cp_index);
			log_B[1] = right_log_likelihood;
			log_q = m_particle.calculate_and_get_adding_binary_log_I_prior_ratio(j, new_trace, trace_index, cp_index);
			
			vector< double > log_Bq(2);
			// create log_Bq's
			log_Bq[0] = log_B[0] + log_q[0];
			log_Bq[1] = log_B[1] + log_q[1];

			double u3 = gsl_ran_flat(r, 0, 1);
			bool accept_cp = log_Bq[1] > log_Bq[0] + log(u3 / (1 - u3));
			if (accept_cp) {
				m_particle.increase_log_likelihood(log_B[1]);
				m_particle.increase_log_binary_I_prior(log_q[1]);
			}
			else {
                m_particle.increase_log_likelihood(log_B[0]);
				m_particle.increase_log_binary_I_prior(log_q[0]);
			}

			if (accept_cp) {
				m_particle.add_new_binary(j, cp_index, right_log_likelihood);
			}
			else {
				m_particle.add_to_binary(j, cp_index, combined_log_likelihood);
			}
		}
	}
	// add a joke binary at the end of the binaries
	m_particle.add_end_binaries(m_number_of_processes);
}

// uses knowledge of the binaries to calculate the log likelihood for the whole particle.
double rj::calculate_total_binary_log_likelihood() {
    double total_log_likelihood = 0;
    if (m_particle.get_dimension() == 0) { // if there are no change points in the model, the likelihood is given by looking at the whole data
        total_log_likelihood += m_pm_ptr->calculate_log_likelihood(m_intercept, m_end + 1);
    }
    for (unsigned int j = 0; j < m_number_of_processes; j++) {
        // calculate the sum of the log likelihood for each segment
        vector< unsigned long int > binary_left_positions = m_particle.get_vector_of_binary_left_positions(j);
        size_t number_of_left_positions = binary_left_positions.size();
        for (unsigned int i = 0; i < number_of_left_positions - 1; i++) {
            double binary_log_likelihood = m_pm_ptr->calculate_log_likelihood(j, binary_left_positions[i], binary_left_positions[i + 1]);
            total_log_likelihood += binary_log_likelihood;
        }
    }
    return total_log_likelihood;
}

void rj::adding_binary_changepoint_setup(const unsigned int & trace_index) {
    unsigned long int lower_position_bound = ((trace_index == 0) ? 0 : m_separators[trace_index - 1]);
    unsigned long int upper_position_bound = ((trace_index == m_number_of_traces - 1) ? m_end + 1 : m_separators[trace_index]);
    m_new_changepoint_position = gsl_rng_uniform_int(r, upper_position_bound - lower_position_bound - 1) + lower_position_bound + 1;
    while (m_particle.does_changepoint_exist_in_particle(m_new_changepoint_position)) {
        m_new_changepoint_position = gsl_rng_uniform_int(r, upper_position_bound - lower_position_bound - 1) + lower_position_bound + 1;
    }
    m_log_proposal_ratio = m_particle.calculate_and_get_add_cp_proposal_ratio(upper_position_bound - lower_position_bound - 1, trace_index, false);
    m_adding_changepoint = changepoint(m_new_changepoint_position);
    unsigned int add_cp_index = m_particle.get_add_changepoint_index();
    // m_log_B contains the log of the Bayes factor L(D1)L(D2) / L(D) for adding the new changepoint to process j
    m_log_B = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_binary_left_log_likelihood = vector< double >(m_number_of_processes, 0);
	m_binary_right_log_likelihood = vector< double >(m_number_of_processes, 0);
    for (unsigned int j = 0; j < m_number_of_processes; j++) {
        // find the index of the cp which begins the binary containing m_tau[cp_index - 1]
        m_binary_left_index = m_particle.get_binary_left_index(j, add_cp_index);
        // calculate the likelihood in the interval [m_tau[m_binary_left_index], new_cp)
        m_binary_left_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_new_changepoint_position);
        // find the index of the cp which begins the binary after the one containing m_tau[cp_index - 1]
        m_binary_right_index = m_particle.get_binary_right_index(j, add_cp_index);
        // calculate the likelihood in the ineterval [new_cp, m_right_index)
        if (m_binary_right_index == m_dimension) {
            m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_new_changepoint_position, m_end + 1);
        } else {
            m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_new_changepoint_position, m_particle.get_changepoint(m_binary_right_index).get_position());
        }
        // calculate the log Bayes factor for adding an effective cp to this process
        m_log_B[j][1] = m_binary_left_log_likelihood[j] + m_binary_right_log_likelihood[j] - m_particle.get_binary_log_likelihood(j, add_cp_index);
    }
    m_log_q = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    for (unsigned int j = 0; j < m_number_of_processes; j++){
        // calculate the log I prior ratio for adding the changepoint and it being ineffective or effective
        m_log_q[j] = m_particle.calculate_and_get_binary_log_I_prior_add_ratio(j);
    }
    // set m_log_Bq = m_log_B + m_log_q
    m_log_Bq = m_log_B;
    for (unsigned int j = 0; j < m_number_of_processes; j++){
        m_log_Bq[j][0] += m_log_q[j][0];
        m_log_Bq[j][1] += m_log_q[j][1];
    }
    // calculate how this move will affect the k prior
    m_log_k_prior_ratio = m_particle.calculate_and_get_add_cp_k_prior_ratio();
    // calculate the proposal ratio for this move (i.e. choosing a birth move and new_position vs choosing a death move and which to kill)
    m_log_acceptance_prob = m_log_k_prior_ratio + m_log_proposal_ratio;
    // calculate log of sum Bq
    for (unsigned int j = 0; j < m_number_of_processes; j++){
        if (m_log_Bq[j][0] > m_log_Bq[j][1]){
            m_log_acceptance_prob += m_log_Bq[j][0] + log(1 + exp(m_log_Bq[j][1] - m_log_Bq[j][0]));
        } else {
            m_log_acceptance_prob += m_log_Bq[j][1] + log(1 + exp(m_log_Bq[j][0] - m_log_Bq[j][1]));
        }
    }
}

void rj::removing_binary_changepoint_setup(const unsigned int & trace_index) {
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    m_log_proposal_ratio = m_particle.calculate_and_get_remove_cp_proposal_ratio(trace_dimension, trace_index, false);
    m_log_B = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_binary_left_log_likelihood = vector< double >(m_number_of_processes, 0);
	m_binary_right_log_likelihood = vector< double >(m_number_of_processes, 0);
    m_binary_merged_log_likelihood = vector< double >(m_number_of_processes, 0);
    for (unsigned int j = 0; j < m_number_of_processes; j++) {
        // find the left index of the binary that contains m_tau[m_h - 1]
        m_binary_left_index = m_particle.get_binary_left_index(j, m_h);
        // find the left index of the binary after the one that contains m_tau[m_h]
        m_binary_right_index = m_particle.get_binary_right_index(j, m_h + 1);
        if (m_particle.get_binary_left_index(j, m_h + 1) == m_h) { // if we are removing an effective changepoint for process j
            // find the log likelihood for the binary containing cp m_(h-1)
            m_binary_left_log_likelihood[j] = m_particle.get_binary_log_likelihood(j, m_h);
            // find the log likelihood for the binary starting at cp m_h
            m_binary_right_log_likelihood[j] = m_particle.get_binary_log_likelihood(j, m_h + 1);
            // calculate the log likelihood of the interval [m_binary_left_index, m_binary_right_index)
            if (m_binary_right_index == m_dimension) {
                m_binary_merged_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_end + 1);
            }
            else {
                m_binary_merged_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_particle.get_changepoint(m_binary_right_index).get_position());
            }
        }
        else { // we are not removing an effective changepoint for process j
            // calculate the log likelihood from this cp to the left of the changepoint m_h to cp m_h
            m_binary_left_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_particle.get_changepoint(m_h).get_position());
            // calculate the log likelihood from m_h to the first effective cp after m_h for process j
            if (m_h == m_dimension - 1){
                m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_h).get_position(), m_end + 1);
            } else {
                if (m_binary_right_index == m_dimension) {
                    m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_h).get_position(), m_end + 1);
                }
                else {
                    m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_h).get_position(), m_particle.get_changepoint(m_binary_right_index).get_position());
                }
            }
            // find the log likelihood of the binary that contains cp m_h - 1
            m_binary_merged_log_likelihood[j] = m_particle.get_binary_log_likelihood(j, m_h);
        }
        m_log_B[j][1] = m_binary_left_log_likelihood[j] + m_binary_right_log_likelihood[j] - m_binary_merged_log_likelihood[j];
    }
    m_log_q = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    for (unsigned int j = 0; j < m_number_of_processes; j++){
        if (m_particle.get_binary_left_index(j, m_h + 1) == m_h) { //if I_{h,j} == 1
            m_log_q[j] = m_particle.calculate_and_get_binary_log_I_prior_remove_ratio(j, true);
        } else {
            m_log_q[j] = m_particle.calculate_and_get_binary_log_I_prior_remove_ratio(j, false);
        }
    }
    // set m_log_Bq = m_log_B + m_log_q
    m_log_Bq = m_log_B;
    for (unsigned int j = 0; j < m_number_of_processes; j++){
        m_log_Bq[j][0] += m_log_q[j][0];
        m_log_Bq[j][1] += m_log_q[j][1];
    }
    m_log_k_prior_ratio = m_particle.calculate_and_get_remove_cp_k_prior_ratio();
    m_log_acceptance_prob = m_log_proposal_ratio + m_log_k_prior_ratio;
    for (unsigned int j = 0; j <  m_number_of_processes; j++){
        if (m_log_Bq[j][0] > m_log_Bq[j][1]){
            m_log_acceptance_prob -= m_log_Bq[j][0] + log(1 + exp(m_log_Bq[j][1] - m_log_Bq[j][0]));
        } else {
            m_log_acceptance_prob -= m_log_Bq[j][1] + log(1 + exp(m_log_Bq[j][0] - m_log_Bq[j][1]));
        }
    }
}

void rj::moving_binary_changepoint_setup(const unsigned int & trace_index) {
    // choose which changepoint to move (we know that there is at least one that can be moved from our definitions of b_k, d_k, etc.)
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    // find the indices of the changepoints before and after m_tau[h] (so that we can sample the position of the new changepoint between them)
	unsigned long int tau_h_minus_1 = m_particle.get_changepoint(m_h - 1).get_position();
	unsigned long int tau_h_plus_1;
	if (m_h == m_dimension - 1){
		tau_h_plus_1 = m_end + 1;
	}
	else {
		tau_h_plus_1 = m_particle.get_changepoint(m_h + 1).get_position();
	}
    // sample the new changepoint position (and don't let m_tau[h] stay in the same position)
	m_new_changepoint_position = gsl_rng_uniform_int(r, tau_h_plus_1 - tau_h_minus_1 - 1) + tau_h_minus_1 + 1;
    if (tau_h_plus_1 - tau_h_minus_1 > 2) {
        while (m_particle.does_changepoint_exist_in_particle(m_new_changepoint_position, m_h, m_h)){
            m_new_changepoint_position = gsl_rng_uniform_int(r, tau_h_plus_1 - tau_h_minus_1 - 1) + tau_h_minus_1 + 1;
        }
    }
	m_log_B = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_log_B_reverse = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_log_q = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_log_q_reverse = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_log_Bq = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_log_Bq_reverse = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
	m_binary_left_log_likelihood = vector< double >(m_number_of_processes, 0);
	m_binary_right_log_likelihood = vector< double >(m_number_of_processes, 0);
    m_binary_merged_log_likelihood = vector< double >(m_number_of_processes, 0);
	m_binary_left_log_likelihood_reverse = vector< double >(m_number_of_processes, 0);
	m_binary_right_log_likelihood_reverse = vector< double >(m_number_of_processes, 0);
	for (unsigned int j = 0; j < m_number_of_processes; j++){
        // find the first effective changepoint to the left of m_tau[m_h]
		m_binary_left_index = m_particle.get_binary_left_index(j, m_h);
        // find the first effective changepoint to the right of m_tau[m_h]
        m_binary_right_index = m_particle.get_binary_right_index(j, m_h + 1);
        // calculate the log likelihood from m_binary_left_index to the new cp position
		m_binary_left_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_new_changepoint_position);
        // calculate the log likelihood from the new cp position to m_binary_right_index
		if (m_binary_right_index == m_dimension){
			m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_new_changepoint_position, m_end + 1);
		}
		else {
			m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_new_changepoint_position, m_particle.get_changepoint(m_binary_right_index).get_position());
		}
        if (m_particle.get_binary_left_index(j, m_h + 1) == m_h) { // if we are moving an effective changepoint
            // need to calculate the merged log likelihood manually as it is not stored in any binary
            if (m_binary_right_index == m_dimension) {
                m_binary_merged_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_end + 1);
            }
            else {
                m_binary_merged_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_particle.get_changepoint(m_binary_right_index).get_position());
            }
        }
        else { // if we are not moving an effeective changepoint
            // merged log likelihood is stored in the binary that containes m_tau[m_h - 1]
            m_binary_merged_log_likelihood[j] = m_particle.get_binary_log_likelihood(j, m_h);
        }
		m_log_B[j][1] = m_binary_left_log_likelihood[j] + m_binary_right_log_likelihood[j] - m_binary_merged_log_likelihood[j];

        if (m_particle.get_binary_left_index(j, m_h + 1) == m_h) { // we are moving an effective changepoint for process j
            // find the log likelihood from m_binary_left_index to the changepoint we are moving
            m_binary_left_log_likelihood_reverse[j] = m_particle.get_binary_log_likelihood(j, m_h);
            // find the log likelihood from m_h to m_binary_right_index
            m_binary_right_log_likelihood_reverse[j] = m_particle.get_binary_log_likelihood(j, m_h + 1);
        }
        else { // we are not removing an effective changepoint for process j
            // calculate the log likelihood from m_binary_left_index to the changepoint we are moving
            m_binary_left_log_likelihood_reverse[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_particle.get_changepoint(m_h).get_position());
            // calculate the log likelihood from the changepoint we are moving to m_binary_right_index
            if (m_binary_right_index == m_dimension) {
                m_binary_right_log_likelihood_reverse[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_h).get_position(), m_end + 1);
            }
            else {
                m_binary_right_log_likelihood_reverse[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_h).get_position(), m_particle.get_changepoint(m_binary_right_index).get_position());
            }
        }
		m_log_B_reverse[j][1] = m_binary_left_log_likelihood_reverse[j] + m_binary_right_log_likelihood_reverse[j] - m_binary_merged_log_likelihood[j];
        
        // calculate m_log_q
		if (m_particle.get_binary_left_index(j, m_h + 1) == m_h) {//if I_{h,j} == 1
			m_log_q[j] = m_log_q_reverse[j] = m_particle.calculate_and_get_binary_log_I_prior_move_ratio(j, true);
		}
		else {
			m_log_q[j] = m_log_q_reverse[j] = m_particle.calculate_and_get_binary_log_I_prior_move_ratio(j, false);
		}
        // m_log_Bq = m_log_B + m_log_q
		m_log_Bq[j][0] = m_log_B[j][0] + m_log_q[j][0];
		m_log_Bq[j][1] = m_log_B[j][1] + m_log_q[j][1];
        // m_log_Bq_reverse = m_log_B_reverse + m_log_q_reverse
		m_log_Bq_reverse[j][0] = m_log_B_reverse[j][0] + m_log_q_reverse[j][0];
		m_log_Bq_reverse[j][1] = m_log_B_reverse[j][1] + m_log_q_reverse[j][1];
	}
    // these are zero because move proposals are symmetric and we are not proposing an increase to the number of changepoints
	m_log_proposal_ratio = 0;
	m_log_k_prior_ratio = 0;
    m_log_acceptance_prob = 0;
    // calculate log of sum of Bq and Bq_reverse
	for (unsigned int j = 0; j < m_number_of_processes; j++){
		if (m_log_Bq[j][0] > m_log_Bq[j][1]){
			m_log_acceptance_prob += m_log_Bq[j][0] + log(1 + exp(m_log_Bq[j][1] - m_log_Bq[j][0]));
		}
		else {
			m_log_acceptance_prob += m_log_Bq[j][1] + log(1 + exp(m_log_Bq[j][0] - m_log_Bq[j][1]));
		}
		if (m_log_Bq_reverse[j][0] > m_log_Bq_reverse[j][1]) {
			m_log_acceptance_prob -= m_log_Bq_reverse[j][0] + log(1 + exp(m_log_Bq_reverse[j][1] - m_log_Bq_reverse[j][0]));
		}
		else {
			m_log_acceptance_prob -= m_log_Bq_reverse[j][1] + log(1 + exp(m_log_Bq_reverse[j][0] - m_log_Bq_reverse[j][1]));
		}
	}
}

void rj::resampling_binary_changepoint_setup(const unsigned int & trace_index) {
    // choose the changepoint index to resample
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    
    m_log_B = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    m_log_B_reverse = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    m_log_q = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    m_log_q_reverse = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    m_log_Bq = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    m_log_Bq_reverse = vector< vector< double > >(m_number_of_processes, vector< double >(2, 0));
    m_binary_left_log_likelihood = vector< double >(m_number_of_processes, 0);
    m_binary_right_log_likelihood = vector< double >(m_number_of_processes, 0);
    m_binary_merged_log_likelihood = vector< double >(m_number_of_processes, 0);
    for (unsigned int j = 0; j < m_number_of_processes; j++){
        // find the first index i to the left of m_h s.t. I_{i,j} = 1
        m_binary_left_index = m_particle.get_binary_left_index(j, m_h);
        //find the first index i to the right of m_h s.t. I_{i,j} = 1
        m_binary_right_index = m_particle.get_binary_right_index(j, m_h + 1);
        if (m_particle.get_binary_left_index(j, m_h + 1) == m_h) { // if we are resampling an effective changepoint for process j
            // find the log likelihood from m_binary_left_index to the changepoint we are moving
            m_binary_left_log_likelihood[j] = m_particle.get_binary_log_likelihood(j, m_h);
            // find the log likelihood from m_h to m_binary_right_index
            m_binary_right_log_likelihood[j] = m_particle.get_binary_log_likelihood(j, m_h + 1);
            // calculate the merged log likelihood
            if (m_binary_right_index == m_dimension) {
                m_binary_merged_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_end + 1);
            }
            else {
                m_binary_merged_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_particle.get_changepoint(m_binary_right_index).get_position());
            }
        }
        else { // if we are resampling an ineffective changepoint for process j
            // calculate the likelihood from the left index to changepint m_h
            m_binary_left_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_particle.get_changepoint(m_h).get_position());
            // calculate the log likelihood in the interval to the right of changepoint m_h
            if (m_binary_right_index == m_dimension){
                m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_h).get_position(), m_end + 1);
            }
            else {
                m_binary_right_log_likelihood[j] = m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_h).get_position(), m_particle.get_changepoint(m_binary_right_index).get_position());
            }
            // calculate the likelihood from m_h to the right_index
            m_binary_merged_log_likelihood[j] = m_particle.get_binary_log_likelihood(j, m_h);
        }
        //recover the likelihood from the left index to the right index (it is given in a binary object)
        m_log_B[j][1] = m_log_B_reverse[j][1] = m_binary_left_log_likelihood[j] + m_binary_right_log_likelihood[j] - m_binary_merged_log_likelihood[j];
        //cout << m_pm_ptr->calculate_log_likelihood(j, m_particle.get_changepoint(m_binary_left_index).get_position(), m_particle.get_changepoint(m_binary_right_index).get_position()) << endl;
        //cout << m_particle.get_binary_left_index(j, m_h + 1) << "\t" << m_h << endl;
        if (m_particle.get_binary_left_index(j, m_h + 1) == m_h){ // if I_{h,j} == 1
            m_log_q[j] = m_log_q_reverse[j] = m_particle.calculate_and_get_binary_log_I_prior_move_ratio(j, true);
        }
        else {
            m_log_q[j] = m_log_q_reverse[j] = m_particle.calculate_and_get_binary_log_I_prior_move_ratio(j, false);
        }
        m_log_Bq[j][0] = m_log_Bq_reverse[j][0] = m_log_B[j][0] + m_log_q[j][0];
        m_log_Bq[j][1] = m_log_Bq_reverse[j][1] = m_log_B[j][1] + m_log_q[j][1];
    }
    // proposal ratio is symmetric
    m_log_proposal_ratio = 0;
    // not adding or removing changepoints
    m_log_k_prior_ratio = 0;
    // move is guaranteed to be accepted
    m_log_acceptance_prob = 0;
}

void rj::binary_acceptance_procedure(const double & u1){
    m_particle.increase_log_k_prior(m_log_k_prior_ratio);
    if (u1 < m_b_k){
        vector< bool > accept_cp;
        for (unsigned int j = 0; j < m_number_of_processes; j++){
            double u3 = gsl_ran_flat(r, 0, 1);
            accept_cp.push_back(m_log_Bq[j][1] > m_log_Bq[j][0] + log(u3 / (1 - u3)));
            if (accept_cp[j]) {
                m_particle.increase_log_likelihood(m_log_B[j][1]);
                m_particle.increase_log_binary_I_prior(m_log_q[j][1]);
            }
            else {
                m_particle.increase_log_binary_I_prior(m_log_q[j][0]);
            }
        }
        m_particle.add_binary_changepoint(m_particle.get_add_changepoint_index(), m_adding_changepoint, accept_cp, m_binary_left_log_likelihood, m_binary_right_log_likelihood);
    }
    else if (u1 < m_d_k){
        vector< bool > remove_effective_cp;
        for (unsigned int j = 0; j < m_number_of_processes; j++){
            if (m_particle.get_binary_left_index(j, m_h + 1) == m_h){//if I_{h,j} == 1
                remove_effective_cp.push_back(true);
                m_particle.increase_log_likelihood(-m_log_B[j][1]);
                m_particle.increase_log_binary_I_prior(-m_log_q[j][1]);
            }
            else {
                remove_effective_cp.push_back(false);
                m_particle.increase_log_binary_I_prior(-m_log_q[j][0]);
            }
        }
        m_particle.remove_binary_changepoint(m_h, remove_effective_cp, m_binary_merged_log_likelihood);
    }
    else if (u1 < m_m_k){
        vector< bool > accept_cp;
        vector< bool > remove_effective_cp;
        for (unsigned int j = 0; j < m_number_of_processes; j++){
            double u3 = gsl_ran_flat(r, 0, 1);
            accept_cp.push_back(m_log_Bq[j][1] > m_log_Bq[j][0] + log(u3 / (1 - u3)));
            remove_effective_cp.push_back(m_particle.get_binary_left_index(j, m_h + 1) == m_h);
            if (accept_cp[j]) {
                if (remove_effective_cp[j]) {
                    m_particle.increase_log_likelihood(m_log_B[j][1] - m_log_B_reverse[j][1]);
                }
                else {
                    m_particle.increase_log_likelihood(m_log_B[j][1]);
                    m_particle.increase_log_binary_I_prior(m_log_q[j][1] - m_log_q_reverse[j][0]);
                }
            }
            else {
                if (remove_effective_cp[j]) {
                    m_particle.increase_log_likelihood(- m_log_B_reverse[j][1]);
                    m_particle.increase_log_binary_I_prior(m_log_q[j][0] - m_log_q_reverse[j][1]);
                }
            }
        }
        m_particle.move_binary_changepoint(m_h, m_new_changepoint_position, remove_effective_cp, accept_cp, m_binary_left_log_likelihood, m_binary_right_log_likelihood, m_binary_merged_log_likelihood);
    }
    else if (u1 < m_r_k){
        vector< bool > accept_cp;
        vector< bool > remove_effective_cp;
        for (unsigned int j = 0; j < m_number_of_processes; j++){
            double u3 = gsl_ran_flat(r, 0, 1);
            accept_cp.push_back(m_log_Bq[j][1] > m_log_Bq[j][0] + log(u3 / (1 - u3)));
            remove_effective_cp.push_back(m_particle.get_binary_left_index(j, m_h + 1) == m_h);
            if (accept_cp[j] && !remove_effective_cp[j]) {
                m_particle.increase_log_likelihood(m_log_B[j][1] - m_log_B_reverse[j][0]);
                m_particle.increase_log_binary_I_prior(m_log_q[j][1] - m_log_q_reverse[j][0]);
            }
            else if (!accept_cp[j] && remove_effective_cp[j]) {
                m_particle.increase_log_likelihood(m_log_B[j][0] - m_log_B_reverse[j][1]);
                m_particle.increase_log_binary_I_prior(m_log_q[j][0] - m_log_q_reverse[j][1]);
            }
        }
        m_particle.resample_binary_changepoint(m_h, remove_effective_cp, accept_cp, m_binary_left_log_likelihood, m_binary_right_log_likelihood, m_binary_merged_log_likelihood);
    }
}

void rj::binary_recording_procedure(){
    m_recorded_binary_dimensions.push_back(m_dimension);
    vector< unsigned long int > changepoint_hist = m_particle.calculate_and_get_binary_changepoint_histogram(m_number_of_changepoint_bins, m_number_of_processes);
    size_t size_of_recorded_changepoints = m_recorded_binary_changepoints.size();
    if (size_of_recorded_changepoints > 0){ //have we started recording changepoint histograms, or is this the first time?
        for (unsigned long int i = 0; i < size_of_recorded_changepoints; i++){
            m_recorded_binary_changepoints[i] += changepoint_hist[i];
        }
    } else {
        m_recorded_binary_changepoints = changepoint_hist;
    }
    double log_posterior = m_particle.get_binary_log_posterior();
    m_recorded_binary_log_posteriors.push_back(log_posterior);
    if (m_binary_MAP_log_posterior < log_posterior){
        m_binary_MAP_particle = m_particle;
        m_binary_MAP_dimension = m_binary_MAP_particle.get_dimension();
        m_binary_MAP_log_posterior = log_posterior;
    }
}

void rj::run_binary_simulation(){
    cout << "starting binary changepoint stage" << endl;
    m_particle.print_likelihood();
    cout << "the posterior is " << m_particle.get_binary_log_posterior() << endl << endl;
    for (unsigned long int iteration = 0; iteration < m_binary_burnin; iteration++) {
        m_dimension = m_particle.get_dimension();
        unsigned int trace_index = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_number_of_traces));
        unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
        if (trace_dimension > 0) {
            m_b_k = 1.0 / 4.0, m_d_k = 2.0 / 4.0, m_m_k = 3.0 / 4.0, m_r_k = 1.0;
        } else {
            m_b_k = 1.0;
        }
        
        //m_particle.check_separator_changepoints();
        //m_particle.print_likelihood();
        //cout << iteration << endl;
        
        double u1 = gsl_ran_flat( r, 0, 1 );
        if (u1 < m_b_k) { //birth
            adding_binary_changepoint_setup(trace_index);
        } else if (u1 < m_d_k) { //death
            removing_binary_changepoint_setup(trace_index);
        } else if (u1 < m_m_k) { //move
            moving_binary_changepoint_setup(trace_index);
        } else if (u1 < m_r_k) { //resample marked vector
            resampling_binary_changepoint_setup(trace_index);
        }
        
        if (m_log_acceptance_prob >= 0 || (log(gsl_ran_flat(r, 0, 1)) < m_log_acceptance_prob)) {
            binary_acceptance_procedure(u1);
        }
        
        //if (m_particle.calculate_and_get_binary_log_posterior(m_number_of_processes) - m_particle.get_binary_log_posterior() > 0.000001 || m_particle.calculate_and_get_binary_log_posterior(m_number_of_processes) - m_particle.get_binary_log_posterior() < -0.000001) {
        //    cout << iteration << '\t' << m_particle.calculate_and_get_binary_log_posterior(m_number_of_processes) - m_particle.get_binary_log_posterior() << endl;
        //}
    }
    
    for (unsigned long int iteration = 0; iteration < m_binary_iterations; iteration++) {
        m_dimension = m_particle.get_dimension();
        unsigned int trace_index = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_number_of_traces));
        unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
        if (trace_dimension > 0) {
            m_b_k = 1.0 / 4.0, m_d_k = 2.0 / 4.0, m_m_k = 3.0 / 4.0, m_r_k = 1.0;
        } else {
            m_b_k = 1.0;
        }
        
        double u1 = gsl_ran_flat(r, 0, 1);
        if (u1 < m_b_k) { //birth
            adding_binary_changepoint_setup(trace_index);
        }
        else if (u1 < m_d_k) { //death
            removing_binary_changepoint_setup(trace_index);
        }
        else if (u1 < m_m_k) { //move
            moving_binary_changepoint_setup(trace_index);
        }
        else if (u1 < m_r_k) { //resample
            resampling_binary_changepoint_setup(trace_index);
        }
        
        if (m_log_acceptance_prob > 0 || (log(gsl_ran_flat(r, 0, 1)) < m_log_acceptance_prob)) {
            binary_acceptance_procedure(u1);
        }
        
        if (m_recording_binary_samples && (iteration % m_binary_thinning == 0)) { //store sample
            binary_recording_procedure();
        }
    }
    cout << "ending binary changepoint stage" << endl;
    m_particle.print_likelihood();
    cout << "the posterior is " << m_particle.get_binary_log_posterior() << endl << endl;
}

void rj::convert_binary_particle_to_full_particle(const double & dirichlet_alpha, const double & rho) {
    m_particle.initiate_regime_vectors(m_number_of_processes);
    /*vector< vector< vector< double > > > sufficient_stats(m_dimension + 2, vector< vector< double > >(m_number_of_processes));
    vector< vector< double > > number_of_observations(m_dimension + 2, vector< double >(m_number_of_processes));
    for (int change = 0; change < m_dimension + 1; change++) {
        m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(change - 1).get_position(), sufficient_stats[change]);
        m_pm_ptr->get_number_of_observations(sufficient_stats[change], number_of_observations[change]);
    }
    m_pm_ptr->get_cumulative_sufficient_data(m_end + 1, sufficient_stats[m_dimension + 1]);
    m_pm_ptr->get_number_of_observations(sufficient_stats[m_dimension + 1], number_of_observations[m_dimension + 1]);*/
    m_particle.set_log_likelihood(0);
    m_particle.set_log_full_I_prior(0);
    m_particle.set_log_regimes_prior(static_cast< double >(m_number_of_processes) * (log(rho) - log(1 - rho)));
    m_particle.set_log_full_separators_prior(0);
    m_particle.set_dirichlet_alpha(dirichlet_alpha);
    m_particle.set_rho(rho);
	set_full_marked_vectors();// m_number_of_processes, sufficient_stats, number_of_observations);
    //check_total_full_log_likelihood(m_particle);
    m_particle.set_log_full_I_prior(m_particle.calculate_and_get_log_full_I_prior(m_number_of_processes));
    m_particle.set_log_regimes_prior(m_particle.calculate_and_get_log_regimes_prior(m_number_of_processes));
    m_particle.set_log_full_separators_prior(m_particle.calculate_and_get_log_full_separators_prior(m_number_of_processes));
    //m_particle.check_full_log_posterior();
}

void rj::set_full_marked_vectors() {//const size_t & number_of_processes), const vector< vector< vector< double > > > & sufficient_statistics, const vector< vector< double > > & number_of_observations) {
    // make sure that the full marked vectors exist for each changepoint so that they can be changed as we go along
    m_particle.set_all_full_marked_vectors_equal_to_binary_marked_vectors(m_number_of_processes);
    for (unsigned int j = 0; j < m_number_of_processes; j++) {
        vector< unsigned long int > binary_left_changepoint_positions = m_particle.get_vector_of_binary_left_positions(j);
        vector< int > binary_left_changepoint_indices = m_particle.get_vector_of_binary_left_indices(j);
        size_t number_of_binaries = m_particle.get_number_of_binaries(j);
        
        unsigned int trace_index = 0;
        // start with binary_index = 0
        unsigned int binary_index = 0;
        unsigned long int left_cp_position = binary_left_changepoint_positions[binary_index];
        //cout << "left cp " << left_cp_position << endl;
        unsigned long int right_cp_position = binary_left_changepoint_positions[binary_index + 1];
        //cout << "right cp " << right_cp_position << endl;
        // we get the sufficient statistics and number of observations for all the processes here when we are only interested in the statistics for process j, but don't want to prematurely optimise here
        vector< vector< double > > stats_1(m_number_of_processes);
        vector< vector< double > > stats_2(m_number_of_processes);
        m_pm_ptr->get_cumulative_sufficient_data(left_cp_position, stats_1);
        m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, stats_2);
        for (size_t i = 0; i < stats_1[j].size(); i++) {
            stats_2[j][i] -= stats_1[j][i];
        }
        double number_of_observations;
        m_pm_ptr->get_number_of_observations(j, stats_2[j], number_of_observations);
        //cout << "number of observations: " << number_of_observations << endl;
        double log_likelihood = m_pm_ptr->calculate_log_likelihood(j, stats_2[j]);
        double log_prior = m_particle.calculate_and_get_full_adding_binary_log_I_prior_ratio(j, 0, binary_left_changepoint_indices[binary_index + 1] - binary_left_changepoint_indices[binary_index], true, trace_index, -1);
        
        vector< unsigned int > transitions = vector< unsigned int >(0);
        // set the regime index to be 0
        unsigned int regime_index = 0;
        vector< unsigned int > transitions_histogram = vector< unsigned int >(1);
        vector< int > right_changepoint_indices;
        for (unsigned int i = binary_left_changepoint_indices[binary_index] + 1; i < binary_left_changepoint_indices[binary_index + 1]; i++) {
            transitions.push_back(regime_index);
            transitions_histogram[regime_index]++;
            right_changepoint_indices.push_back(i);
        }
        transitions.push_back(-1);
        right_changepoint_indices.push_back(binary_left_changepoint_indices[binary_index + 1]);
        // add this regime to the particle
        m_particle.add_new_regime(j, right_changepoint_indices, transitions, transitions_histogram, stats_2[j], number_of_observations, log_likelihood, m_number_of_traces, true, trace_index, -1);
        m_particle.increase_log_likelihood(log_likelihood);
        m_particle.increase_log_full_I_prior(log_prior, true, j, trace_index, true);
        //m_particle.check_transitions_out();
        //cout << m_particle.calculate_and_get_log_full_I_prior(j+1) - m_particle.get_log_full_I_prior() << endl;
        
        for (unsigned int binary_index = 1; binary_index < number_of_binaries - 1; binary_index++) {
            unsigned long int left_cp_position = binary_left_changepoint_positions[binary_index];
            unsigned long int right_cp_position = binary_left_changepoint_positions[binary_index + 1];
            int left_cp_index = binary_left_changepoint_indices[binary_index];
            int right_cp_index = binary_left_changepoint_indices[binary_index + 1];
            bool new_trace = false;
            if (m_particle.is_changepoint_index_separator_index(left_cp_index)) {
                trace_index++;
                new_trace = true;
            }
            unsigned int previous_regime = m_particle.get_previous_regime(left_cp_index, j);
            
            // we get the sufficient statistics and number of observations for all the processes here when we are only interested in the statistics for process j, but don't want to prematurely optimise here
            vector< vector< double > > stats_1(m_number_of_processes);
            vector< vector< double > > stats_2(m_number_of_processes);
            m_pm_ptr->get_cumulative_sufficient_data(left_cp_position, stats_1);
            m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, stats_2);
            for (size_t i = 0; i < stats_1[j].size(); i++) {
                stats_2[j][i] -= stats_1[j][i];
            }
            double number_of_observations;
            m_pm_ptr->get_number_of_observations(j, stats_2[j], number_of_observations);
            double log_likelihood = m_pm_ptr->calculate_log_likelihood(j, stats_2[j]);
            
            size_t number_of_regimes = m_particle.get_number_of_regimes(j);
            vector< double > log_B(number_of_regimes + 1, 0);
            vector< double > log_likelihoods_with_right_sufficient_statistics(number_of_regimes + 1, 0);
            vector< double > log_q(number_of_regimes + 1, 0);
            // now choose the regime to which we will add this binary
            for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
                //cout << "reg suff stat: " << regime_sufficient_stats[1] << endl;
                // calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
                vector< double > regime_sufficient_stats_with_right_sufficient_statistics = regime_sufficient_stats;
                for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
                    regime_sufficient_stats_with_right_sufficient_statistics[index] += stats_2[j][index];
                }
                double regime_log_likelihood = m_particle.get_regime_log_likelihood(j, regime);
                //cout << "regime log likelihood: " << regime_log_likelihood << endl;
                double regime_log_likelihood_with_sufficient_stats = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics);
                //cout << "regime log likelihood with sufficient stats: " << regime_log_likelihood_with_sufficient_stats << endl;
                log_B[regime] = regime_log_likelihood_with_sufficient_stats - regime_log_likelihood;
                log_likelihoods_with_right_sufficient_statistics[regime] = regime_log_likelihood_with_sufficient_stats;
                log_q[regime] = m_particle.calculate_and_get_full_adding_binary_log_I_prior_ratio(j, regime, right_cp_index - left_cp_index, new_trace, trace_index, previous_regime);
            }
            // also propose adding a new regime
            double regime_log_likelihood_with_sufficient_stats = log_likelihood;
            log_B[number_of_regimes] = regime_log_likelihood_with_sufficient_stats;
            log_likelihoods_with_right_sufficient_statistics[number_of_regimes] = regime_log_likelihood_with_sufficient_stats;
            log_q[number_of_regimes] = m_particle.calculate_and_get_full_adding_binary_log_I_prior_ratio(j, static_cast< unsigned int >(number_of_regimes), right_cp_index - left_cp_index, new_trace, trace_index, previous_regime);
            
            vector< double > log_Bq(number_of_regimes + 1);
            // create log_Bq's
            for (unsigned int regime = 0; regime < number_of_regimes + 1; regime++) {
                log_Bq[regime] = log_B[regime] + log_q[regime];
            }
            // calculate the ordering of the elements of m_log_Bq so that we can use fancy log and exp tricks to calculate the sum of the log Bq values.
            vector< unsigned int > log_Bq_descending_order(number_of_regimes + 1, 0);
            for (unsigned int i = 1; i < number_of_regimes + 1; i++) {
                log_Bq_descending_order[i] = i;
            }
            calculate_vector_descending_order(number_of_regimes + 1, log_Bq, log_Bq_descending_order);
            
            double log_of_sum_Bq = log_Bq[log_Bq_descending_order[0]];
            double temp_sum = 1;
            double largest_Bq = log_Bq[log_Bq_descending_order[0]];
            for (unsigned int index = 1; index < number_of_regimes + 1; index++){
                temp_sum += exp(log_Bq[log_Bq_descending_order[index]] - largest_Bq);
            }
            log_of_sum_Bq += log(temp_sum);
            
            // choose to which regime we will add this binary
            temp_sum = 0;
            double regime_chooser = log(gsl_ran_flat(r, 0, 1)) + log_of_sum_Bq;
            unsigned int index = 0;
            bool chosen = false;
            do {
                temp_sum += exp(log_Bq[log_Bq_descending_order[index]] - log_Bq[log_Bq_descending_order[0]]);
                chosen = regime_chooser <= log_Bq[log_Bq_descending_order[0]] + log(temp_sum);
                index++;
            }
            while(!chosen);
            unsigned int new_regime = log_Bq_descending_order[index - 1];
            m_particle.increase_log_likelihood(log_B[new_regime]);
            //cout << "log_B[new_regime]: " << log_B[new_regime] << endl;
            bool adding_new_regime = new_regime == number_of_regimes;
            m_particle.increase_log_full_I_prior(log_q[new_regime], adding_new_regime, j, trace_index, new_trace); // if a new regime is added, the (1-rho) in m_log_q[j][new_regimes[j]] will be added to m_log_regime_prior and taken away from m_log_full_I_prior
            
            vector< unsigned int > transitions = vector< unsigned int >(0);
            // set the regime index to be the new regime number
            unsigned int regime_index = new_regime;
            vector< unsigned int > transitions_histogram = vector< unsigned int >(number_of_regimes + (adding_new_regime ? 1 : 0));
            vector< int > right_changepoint_indices(0);
            for (int i = left_cp_index + 1; i < right_cp_index; i++) {
                transitions.push_back(regime_index);
                transitions_histogram[regime_index]++;
                right_changepoint_indices.push_back(i);
            }
            transitions.push_back(-1);
            right_changepoint_indices.push_back(right_cp_index);
            
            if (adding_new_regime) {
                m_particle.add_new_regime(j, right_changepoint_indices, transitions, transitions_histogram, stats_2[j], number_of_observations, log_likelihood, m_number_of_traces, new_trace, trace_index, previous_regime);
            } else {
                m_particle.add_binary_to_regime(j, regime_index, right_changepoint_indices, transitions, transitions_histogram, stats_2[j], number_of_observations, log_likelihoods_with_right_sufficient_statistics[regime_index], m_number_of_traces, new_trace, trace_index, previous_regime);
            }
            //m_particle.check_transitions_out();
            //cout << m_particle.calculate_and_get_log_full_I_prior(j+1) - m_particle.get_log_full_I_prior() << endl;
        }
    }
    // set the number of unobserved regimes for each process to be 0
    m_particle.set_all_regimes_to_be_observed(m_number_of_processes);
}

// calculate the full log likelihood from the changepoint positions
void rj::check_total_full_log_likelihood(particle & P) {
    double log_likelihood = 0;
    for (unsigned int process = 0; process < m_number_of_processes; process++) {
        size_t number_of_regimes = P.get_number_of_regimes(process);
        for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
            double regime_log_likelihood = 0;
            vector< int > regime_right_changepoint_indices = P.get_right_changepoint_indices(process, regime);
            vector< vector< double > > regime_sufficient_stats(m_number_of_processes);
            if (0 < regime_right_changepoint_indices.size()) {
                unsigned long int left_cp_position = P.get_changepoint(regime_right_changepoint_indices[0] - 1).get_position();
                unsigned long int right_cp_position;
                if (regime_right_changepoint_indices[0] == P.get_dimension()) {
                    right_cp_position = m_end + 1;
                }
                else {
                    right_cp_position = P.get_changepoint(regime_right_changepoint_indices[0]).get_position();
                }
                vector< vector< double > > stats_1(m_number_of_processes);
                m_pm_ptr->get_cumulative_sufficient_data(left_cp_position, stats_1);
                m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, regime_sufficient_stats);
                for (size_t i = 0; i < stats_1[process].size(); i++) {
                    regime_sufficient_stats[process][i] -= stats_1[process][i];
                }
            }
            for (unsigned int index = 1; index < regime_right_changepoint_indices.size(); index++) {
                // get sufficient stats for the interval from CP regime_right_changepoint_indices[index] - 1 to CP regime_right_changepoint_indices[index]
                unsigned long int left_cp_position = P.get_changepoint(regime_right_changepoint_indices[index] - 1).get_position();
                unsigned long int right_cp_position;
                if (regime_right_changepoint_indices[index] == P.get_dimension()) {
                    right_cp_position = m_end + 1;
                }
                else {
                    right_cp_position = P.get_changepoint(regime_right_changepoint_indices[index]).get_position();
                }
                vector< vector< double > > stats_1(m_number_of_processes);
                vector< vector< double > > stats_2(m_number_of_processes);
                m_pm_ptr->get_cumulative_sufficient_data(left_cp_position, stats_1);
                m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, stats_2);
                for (size_t i = 0; i < stats_1[process].size(); i++) {
                    stats_2[process][i] -= stats_1[process][i];
                    regime_sufficient_stats[process][i] += stats_2[process][i];
                }
            }
            if (0 < regime_right_changepoint_indices.size()) {
                regime_log_likelihood = m_pm_ptr->calculate_log_likelihood(process, regime_sufficient_stats[process]);
            } // else leave it equal to 0 because this regime is unobserved
            // check that the regime_log_likelihood we have calculated equals the log likelihood stored in the regime
            if (abs(regime_log_likelihood - P.get_regime_log_likelihood(process, regime)) > 0.00001) {
                cout << "log likelihood doesn't match" << endl;
            }
            log_likelihood += regime_log_likelihood;
        }
    }
    if (abs(log_likelihood - P.get_log_likelihood()) > 0.000001) {
        cout << "total log likelihood doesn't match" << endl;
    }
}

/*void rj::check_adding_changepoint(const unsigned int & add_cp_index) {
    // cycing through different regimes to add to the particle, check if the log likelihood ratio is equal to what m_log_B implies it would be and m_log_q values are right as well
    // work out which process has the most regimes
    
    for (unsigned int process = 0; process < m_number_of_processes; process++) {
        size_t number_of_regimes = m_particle.get_number_of_regimes(process);
        for (unsigned int regime = 0; regime < number_of_regimes + 1; regime++) {
            // create a copy of the particle
            particle P = m_particle;
            // for log likelihood need to have the changepoint added and the right_changepoint indices for the regime we are adding to and the previous one need to be amended
            changepoint Frederick(m_new_changepoint_position);
            P.add_full_changepoint(add_cp_index, Frederick, <#const vector<unsigned int> &new_regimes#>, <#const vector<vector<double> > &sufficient_statistics#>, <#const vector<vector<double> > &log_likelihoods#>, <#const vector<double> &previous_log_likelihoods#>, <#const vector<double> &number_of_observations#>)
            if (m_particle.is_changepoint_index_separator_index(add_cp_index - 1)) { // if the previous
                
            }
            // for log_full_I_prior need to alter the right transitions histogram for the regime we are adding to and the previous regime.
            
        }
    }
}*/

// A and order both have size length. order contains 0, 1, ..., length - 1, and A is a copy of a (usually) unordered vector. order will be set to be the descending ordering of A, so order[0] will give the index of A that is largest. The ordering will not detect ties and will simply order whichever is first in the vector as the largest if there is a tie.
void rj::calculate_vector_descending_order(const size_t & length, vector< double >  A, vector< unsigned int > & order) {
    int i, j;
    unsigned int b;
    double a;
    for (j = 1; j < length; j++) {
        a = A[j];
        b = order[j];
        i = j - 1;
        while (i >= 0 && A[i] < a) {
            A[i + 1] = A[i];
            order[i + 1] = order[i];
            i--;
        }
        A[i + 1] = a;
        order[i + 1] = b;
    }
}

void rj::adding_full_changepoint_setup(const unsigned int & trace_index) {
    unsigned long int lower_position_bound = ((trace_index == 0) ? 0 : m_separators[trace_index - 1]);
    unsigned long int upper_position_bound = ((trace_index == m_number_of_traces - 1) ? m_end + 1 : m_separators[trace_index]);
    m_new_changepoint_position = gsl_rng_uniform_int(r, upper_position_bound - lower_position_bound - 1) + lower_position_bound + 1;
    while (m_particle.does_changepoint_exist_in_particle(m_new_changepoint_position)) {
        m_new_changepoint_position = gsl_rng_uniform_int(r, upper_position_bound - lower_position_bound - 1) + lower_position_bound + 1;
    }
    m_log_proposal_ratio = m_particle.calculate_and_get_add_cp_proposal_ratio(upper_position_bound - lower_position_bound - 1, trace_index, false);
    m_adding_changepoint = changepoint(m_new_changepoint_position);
    unsigned int add_cp_index = m_particle.get_add_changepoint_index();
    // calculate how this move will affect the k prior
    m_log_k_prior_ratio = m_particle.calculate_and_get_add_cp_k_prior_ratio();
    m_log_acceptance_prob = m_log_k_prior_ratio + m_log_proposal_ratio;
	unsigned long int right_cp_position;
	if (add_cp_index == m_dimension) {
		right_cp_position = m_end + 1;
	}
	else {
		right_cp_position = m_particle.get_changepoint(add_cp_index).get_position();
	}
	unsigned long int left_cp_position = m_particle.get_changepoint(add_cp_index - 1).get_position();
	// m_log_B contains the log of the Bayes factors for adding a changepoint with each regime, so m_log_B[i] = log_Bayes_factor_if_I_h_j_prime_equals_i
	m_log_B = vector< vector< double > >(0);
	// m_log_q contains the log of the I prior ratio for adding a changepoint with each regime, so m_log_B[i] = log_I_prior_ratio_if_I_h_j_prime_equals_i
	m_log_q = vector< vector< double > >(0);
    // m_log_Bq contains the log of (the I prior ratio times the likelihood)
    m_log_Bq = vector< vector< double > >(0);
    // m_log_Bq_ordering gives the order of the elements of m_log_Bq. So m_log_Bq_ordering[j][0] gives the index of the largest element in m_log_Bq[j], etc
    m_log_Bq_descending_order = vector< vector< unsigned int > >(0);
    // m_log_of_sum_Bq contains the log of the sum of the Bq values for each process (useful for calculating which regime to choose and for acceptance probability calculations)
    m_log_of_sum_Bq = vector< double >(m_number_of_processes, 0);
    //m_previous_log_likelihoods_without_right_sufficient_statistics gives (for each process) the likelihood for the previous regime without the right sufficient stats
    m_previous_log_likelihoods_without_right_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    // m_log_likelihoods_with_right_sufficient_statistics gives (for each process) the likelihood for each regime if the right sufficient stats are added to it. Equals 0 for the previous regime for each process, as the 'adding a changepoint' procedure knows that it shouldn't be setting new log likelihoods if new regime == previous regime
    m_log_likelihoods_with_right_sufficient_statistics = vector< vector< double > >(0);
	// m_left_sufficient_statistics holds the sufficient statistics for the interval from (the changepoint prior to m_new_changepoint_position) to m_new_changepoint_position for process j.
	// the assignment here looks wrong, but it will have the cumulative data up to left_cp_position subtracted later
    m_left_sufficient_statistics = vector< vector< double > >(m_number_of_processes);
	m_pm_ptr->get_cumulative_sufficient_data(m_new_changepoint_position, m_left_sufficient_statistics);
	// m_right_sufficient_statistics holds the sufficient statistics for the interval from m_new_changepoint_position to (the changepoint after m_new_changepoint_position) for process j
	// the assignment here looks wrong, but it will have the cumulative data up to m_new_cp_position subtracted later
    m_right_sufficient_statistics = vector< vector< double > >(m_number_of_processes);
	m_pm_ptr->get_cumulative_sufficient_data(right_cp_position, m_right_sufficient_statistics);
    vector< vector< double > > sufficient_statistics_up_to_left_cp_position(m_number_of_processes);
    m_pm_ptr->get_cumulative_sufficient_data(left_cp_position, sufficient_statistics_up_to_left_cp_position);
	for (unsigned int j = 0; j < m_number_of_processes; j++) {
		size_t number_of_regimes = m_particle.get_number_of_regimes(j);
		m_log_B.push_back(vector< double >(number_of_regimes + 1, 0));
		m_log_q.push_back(vector< double >(number_of_regimes + 1, 0));
        m_log_likelihoods_with_right_sufficient_statistics.push_back(vector< double >(number_of_regimes + 1, 0));
		for (unsigned int index = 0; index < m_right_sufficient_statistics[j].size(); index++) {
			// this is correct, as m_left_sufficient_statistics currently holds the information up to the new changepoint position
			m_right_sufficient_statistics[j][index] -= m_left_sufficient_statistics[j][index];
			// we now make m_left_sufficient_statistics correct
			m_left_sufficient_statistics[j][index] -= sufficient_statistics_up_to_left_cp_position[j][index];
		}
		// get the sufficient statistics for the regime that affects the changepoint prior to the new changepoint
		unsigned int previous_regime = m_particle.get_previous_regime(add_cp_index, j);
		vector< double > previous_regime_sufficient_stats = m_particle.get_sufficient_statistics(j, previous_regime);
		// calculate the likelihood for the regime that affects the changepoint prior to the new changepoint
		double previous_regime_log_likelihood = m_particle.get_regime_log_likelihood(j, previous_regime);
		// calculate the likelihood for the previous regime with the right sufficient statistics removed
		vector< double > previous_regime_sufficient_statistics_without_right_sufficient_statistics = previous_regime_sufficient_stats;
		for (unsigned int index = 0; index < previous_regime_sufficient_statistics_without_right_sufficient_statistics.size(); index++) {
			previous_regime_sufficient_statistics_without_right_sufficient_statistics[index] -= m_right_sufficient_statistics[j][index];
		}
		double previous_regime_log_likelihood_without_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_right_sufficient_statistics);
        m_previous_log_likelihoods_without_right_sufficient_statistics[j] = previous_regime_log_likelihood_without_right_sufficient_statistics;
        // for each regime calculate the Bayes factor for assigning this regime to the new changepoint.
		for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
			if (regime != previous_regime) { // check if there is anything to calculate. If the regime we are considering for the new changepoint is the same as the previous one, then the Bayes factor is 1.
				vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
				// calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
				vector< double > regime_sufficient_stats_with_right_sufficient_statistics = m_particle.get_sufficient_statistics(j, regime);
				for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
					regime_sufficient_stats_with_right_sufficient_statistics[index] += m_right_sufficient_statistics[j][index];
				}
                double regime_log_likelihood_with_sufficient_stats = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics);
				m_log_B[j][regime] = previous_regime_log_likelihood_without_right_sufficient_statistics + regime_log_likelihood_with_sufficient_stats - previous_regime_log_likelihood - m_particle.get_regime_log_likelihood(j, regime);
                m_log_likelihoods_with_right_sufficient_statistics[j][regime] = regime_log_likelihood_with_sufficient_stats;
			}
		}
		// also propose adding a new regime
        double regime_log_likelihood_with_sufficient_stats = m_pm_ptr->calculate_log_likelihood(j, m_right_sufficient_statistics[j]);
        m_log_B[j][number_of_regimes] = previous_regime_log_likelihood_without_right_sufficient_statistics + regime_log_likelihood_with_sufficient_stats - previous_regime_log_likelihood;
        m_log_likelihoods_with_right_sufficient_statistics[j][number_of_regimes] = regime_log_likelihood_with_sufficient_stats;
        
        // for each regime calculate marked vector prior ratio
        if (add_cp_index == m_dimension || m_particle.is_changepoint_index_separator_index(add_cp_index)) {
            for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                m_log_q[j][regime] = m_particle.full_log_I_prior_add_ratio(add_cp_index, false, j, previous_regime, regime);
            }
            m_log_q[j][number_of_regimes] = m_particle.full_log_I_prior_add_ratio(add_cp_index, true, j, previous_regime, static_cast< unsigned int >(number_of_regimes));
        }
        else {
            unsigned int subsequent_regime = m_particle.get_previous_regime(add_cp_index + 1, j);
            for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                m_log_q[j][regime] = m_particle.full_log_I_prior_add_ratio(add_cp_index, false, j, previous_regime, regime, subsequent_regime);
            }
            m_log_q[j][number_of_regimes] = m_particle.full_log_I_prior_add_ratio(add_cp_index, true, j, previous_regime, static_cast< unsigned int >(number_of_regimes), subsequent_regime);
        }
        
        //check_adding_changepoint(add_cp_index);
        
        // set m_log_Bq = m_log_B + m_log_q
        m_log_Bq.push_back(vector< double >(number_of_regimes + 1, 0));
        for (unsigned int regime = 0; regime < number_of_regimes + 1; regime++){
            m_log_Bq[j][regime] = m_log_B[j][regime] + m_log_q[j][regime];
        }
        
        // check if regime r_j is unobserved: if so, the Bq term will equal 0 as there would be no reverse move.
        if (m_particle.is_regime_unobserved(j, number_of_regimes - 1)) {
            m_log_Bq[j][number_of_regimes - 1] = -1e300;
        }
        
        // calculate the ordering of the elements of m_log_Bq so that we can use fancy log and exp tricks to calculate the sum of the log Bq values.
        m_log_Bq_descending_order.push_back(vector< unsigned int >(number_of_regimes + 1, 0));
        for (unsigned int i = 1; i < number_of_regimes + 1; i++) {
            m_log_Bq_descending_order[j][i] = i;
        }
        calculate_vector_descending_order(number_of_regimes + 1, m_log_Bq[j], m_log_Bq_descending_order[j]);
        // calculate log of sum Bq
        double largest_Bq = m_log_Bq[j][m_log_Bq_descending_order[j][0]];
        m_log_of_sum_Bq[j] = m_log_Bq[j][m_log_Bq_descending_order[j][0]];
        double temp_sum = 1;
        for (unsigned int index = 1; index < number_of_regimes + 1; index++){
            temp_sum += exp(m_log_Bq[j][m_log_Bq_descending_order[j][index]] - largest_Bq);
        }
        m_log_of_sum_Bq[j] += log(temp_sum);
        m_log_acceptance_prob += m_log_of_sum_Bq[j];
    }
}

void rj::removing_full_changepoint_setup(const unsigned int & trace_index) {
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    m_log_proposal_ratio = m_particle.calculate_and_get_remove_cp_proposal_ratio(trace_dimension, trace_index, false);
    m_log_k_prior_ratio = m_particle.calculate_and_get_remove_cp_k_prior_ratio();
    m_log_acceptance_prob = m_log_k_prior_ratio + m_log_proposal_ratio;
    
    // m_log_B_reverse[i] contains log Bayes factor if I_h_j_prime = i vs no changepoint here. m_log_q_reverse is the log I ratio for the same setting. m_log_Bq_reverse is their sum
    m_log_B_reverse = vector< vector< double > >(0);
    m_log_q_reverse = vector< vector< double > >(0);
    m_log_Bq_reverse = vector< vector< double > >(0);
    // m_log_Bq_reverse_descending_order[j][0] gives the index of the largest element in m_log_Bq_reverse, m_log_of_sum_Bq_reverse[j] gives the log of the sum of the Bq_reverse[j]
    m_log_Bq_reverse_descending_order = vector< vector< unsigned int > >(0);
    m_log_of_sum_Bq_reverse = vector< double >(m_number_of_processes, 0);
    //m_previous_log_likelihoods_with_right_sufficient_statistics_reverse gives (for each process) the likelihood for the previous regime with the reverse right sufficient stats
    m_previous_log_likelihoods_with_right_sufficient_statistics_reverse = vector< double >(m_number_of_processes, 0);
    // m_actual_log_likelihoods_without_right_sufficient_statistics_reverse gives (for each process) the likelihood for the removed regime if the right sufficient stats are removed from it. Equals 0 for the previous regime for each process, as the 'removing a changepoint' procedure knows that it shouldn't be setting new log likelihoods if new regime == previous regime
    m_actual_log_likelihoods_without_right_sufficient_statistics_reverse = vector< double >(m_number_of_processes, 0);
    // m_right_sufficient_statistics_reverse holds the sufficient statistics for the interval from m_h to m_h+1 for process j
    // the assignment here looks wrong, but it will have the cumulative data up to m_h subtracted later
    m_right_sufficient_statistics_reverse = vector< vector< double > >(m_number_of_processes);
    // calculate how this move will affect the k prior
    if (m_h == m_dimension - 1) {
        m_pm_ptr->get_cumulative_sufficient_data(m_end + 1, m_right_sufficient_statistics_reverse);
    }
    else {
        m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h + 1).get_position(), m_right_sufficient_statistics_reverse);
    }
    m_left_sufficient_statistics_reverse = vector< vector< double > >(m_number_of_processes); // Is this needed?
    m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h).get_position(), m_left_sufficient_statistics_reverse);
    vector< vector< double > > sufficient_statistics_up_to_left_cp_position(m_number_of_processes);
    m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h - 1).get_position(), sufficient_statistics_up_to_left_cp_position);
    m_removing_unobserved_regimes = vector< bool >(m_number_of_processes, false);
    for (unsigned int j = 0; j < m_number_of_processes; j++) {
        size_t number_of_regimes = m_particle.get_number_of_regimes(j);
        // get the actual regime for changepoint m_h
        unsigned int actual_regime = m_particle.get_previous_regime(m_h + 1, j);
        // if we are removing a whole regime from the process then we need to subtract 1 from the number of regimes.
        // if deleting tau_h means that regime r_j becomes unobserved then the number of regimes is reduced by 1.
        if (m_particle.removing_full_changepoint_leaves_highest_regime_unobserved(j, actual_regime)) {
            number_of_regimes--;
            m_removing_unobserved_regimes[j] = true;
        }
        m_log_B_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
        m_log_q_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
        for (unsigned int index = 0; index < m_right_sufficient_statistics_reverse[j].size(); index++) {
            // this is correct, as m_left_sufficient_statistics_reverse currently holds the information up to m_h
            m_right_sufficient_statistics_reverse[j][index] -= m_left_sufficient_statistics_reverse[j][index];
            // we now make m_left_sufficient_statistics_reverse correct
            m_left_sufficient_statistics_reverse[j][index] -= sufficient_statistics_up_to_left_cp_position[j][index];
        }
        // get the sufficient statistics for the regime that affects the changepoint prior to the new changepoint
        unsigned int previous_regime = m_particle.get_previous_regime(m_h, j);
        vector< double > previous_regime_sufficient_stats = m_particle.get_sufficient_statistics(j, previous_regime);
        // calculate the likelihood for the regime that affects the changepoint prior to the new changepoint
        double previous_regime_log_likelihood = m_particle.get_regime_log_likelihood(j, previous_regime);
        // calculate the likelihood for the previous regime with the right sufficient statistics removed
        vector< double > previous_regime_sufficient_statistics_with_right_sufficient_statistics = previous_regime_sufficient_stats;
        for (unsigned int index = 0; index < previous_regime_sufficient_statistics_with_right_sufficient_statistics.size(); index++) {
            previous_regime_sufficient_statistics_with_right_sufficient_statistics[index] += m_right_sufficient_statistics_reverse[j][index];
        }
        double previous_regime_log_likelihood_with_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_with_right_sufficient_statistics);
        m_previous_log_likelihoods_with_right_sufficient_statistics_reverse[j] = previous_regime_log_likelihood_with_right_sufficient_statistics;
        // for each regime calculate the reverse Bayes factor for assigning this regime to the new changepoint.
        for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
            if (previous_regime == actual_regime) {
                if (regime == previous_regime) {
                    m_log_B_reverse[j][regime] = 0;
                }
                else {
                    vector< double > previous_regime_sufficient_statistics_without_right_sufficient_statistics = previous_regime_sufficient_stats;
                    for (unsigned int index = 0; index < previous_regime_sufficient_statistics_with_right_sufficient_statistics.size(); index++) {
                        previous_regime_sufficient_statistics_without_right_sufficient_statistics[index] -= m_right_sufficient_statistics_reverse[j][index];
                    }
                    double previous_regime_log_likelihood_without_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_right_sufficient_statistics);
                    vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
                    // calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
                    vector< double > regime_sufficient_stats_with_right_sufficient_statistics = m_particle.get_sufficient_statistics(j, regime);
                    for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
                        regime_sufficient_stats_with_right_sufficient_statistics[index] += m_right_sufficient_statistics_reverse[j][index];
                    }
                    m_log_B_reverse[j][regime] = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics) + previous_regime_log_likelihood_without_right_sufficient_statistics - m_particle.get_regime_log_likelihood(j, regime) - previous_regime_log_likelihood;
                }
            }
            else {
                if (regime == previous_regime) {
                    m_log_B_reverse[j][regime] = 0;
                }
                else {
                    if (regime == actual_regime) {
                        vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
                        // calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
                        vector< double > regime_sufficient_stats_without_right_sufficient_statistics = m_particle.get_sufficient_statistics(j, regime);
                        for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
                            regime_sufficient_stats_without_right_sufficient_statistics[index] -= m_right_sufficient_statistics_reverse[j][index];
                        }
                        double regime_log_likelihood_without_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_without_right_sufficient_statistics);
                        m_actual_log_likelihoods_without_right_sufficient_statistics_reverse[j] = regime_log_likelihood_without_right_sufficient_statistics;
                        m_log_B_reverse[j][regime] = m_particle.get_regime_log_likelihood(j, regime) + previous_regime_log_likelihood - regime_log_likelihood_without_right_sufficient_statistics - previous_regime_log_likelihood_with_right_sufficient_statistics;
                    }
                    else {
                        vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
                        // calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
                        vector< double > regime_sufficient_stats_with_right_sufficient_statistics = m_particle.get_sufficient_statistics(j, regime);
                        for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
                            regime_sufficient_stats_with_right_sufficient_statistics[index] += m_right_sufficient_statistics_reverse[j][index];
                        }
                        m_log_B_reverse[j][regime] = m_pm_ptr->calculate_log_likelihood(j,regime_sufficient_stats_with_right_sufficient_statistics) + previous_regime_log_likelihood - previous_regime_log_likelihood_with_right_sufficient_statistics - m_particle.get_regime_log_likelihood(j, regime);
                    }
                }
            }
        }
        // also propose adding a new regime
        if (previous_regime == actual_regime) {
            // calculate the likelihood for the previous regime with right interval removed
            vector< double > previous_regime_sufficient_statistics_without_right_sufficient_statistics = previous_regime_sufficient_stats;
            for (unsigned int index = 0; index < previous_regime_sufficient_statistics_with_right_sufficient_statistics.size(); index++) {
                previous_regime_sufficient_statistics_without_right_sufficient_statistics[index] -= m_right_sufficient_statistics_reverse[j][index];
            }
            double previous_regime_log_likelihood_without_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_right_sufficient_statistics);
            m_log_B_reverse[j][number_of_regimes] = previous_regime_log_likelihood_without_right_sufficient_statistics + m_pm_ptr->calculate_log_likelihood(j, m_right_sufficient_statistics_reverse[j]) - previous_regime_log_likelihood;
        }
        else {
            m_log_B_reverse[j][number_of_regimes] = previous_regime_log_likelihood + m_pm_ptr->calculate_log_likelihood(j, m_right_sufficient_statistics_reverse[j]) - previous_regime_log_likelihood_with_right_sufficient_statistics;
        }
        
        // for each regime calculate marked vector prior ratio
        if (m_h == m_dimension - 1 || m_particle.is_changepoint_index_separator_index(m_h + 1)) {
            for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                m_log_q_reverse[j][regime] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, regime, false, actual_regime, m_removing_unobserved_regimes[j]);
            }
            m_log_q_reverse[j][number_of_regimes] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, static_cast< unsigned int >(number_of_regimes), true, actual_regime, m_removing_unobserved_regimes[j]);
        }
        else {
            unsigned int subsequent_regime = m_particle.get_previous_regime(m_h + 2, j);
            for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                m_log_q_reverse[j][regime] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, regime, false, actual_regime, m_removing_unobserved_regimes[j], subsequent_regime);
            }
            m_log_q_reverse[j][number_of_regimes] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, static_cast< unsigned int >(number_of_regimes), true, actual_regime, m_removing_unobserved_regimes[j], subsequent_regime);
        }
        // calculate the trace index
        m_particle.is_changepoint_index_separator_index(m_h); // only running this to set the trace_index. If not run, m_trace_index may give the trace index of the next trace because is_changepoint_index_separator_index(index + 1) has just been run
        
        // set m_log_Bq = m_log_B + m_log_q
        m_log_Bq_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
        for (unsigned int regime = 0; regime < number_of_regimes + 1; regime++){
            m_log_Bq_reverse[j][regime] = m_log_B_reverse[j][regime] + m_log_q_reverse[j][regime];
        }
        // check if regime number_of_regimes - 1 is unobserved. If removing tau_h reduces the number of regimes by 1, this will already be accounted for
        if (m_particle.is_regime_unobserved(j, number_of_regimes - 1)) {
            m_log_Bq_reverse[j][number_of_regimes - 1] = -1e300;
        }
        // initialise m_log_Bq_reverse_descending_order and then calculate the ordering of the elements of m_log_Bq so that we can use fancy log and exp tricks to calculate the sum of the log Bq values.
        m_log_Bq_reverse_descending_order.push_back(vector< unsigned int >(number_of_regimes + 1, 0));
        for (unsigned int i = 1; i < number_of_regimes + 1; i++) {
            m_log_Bq_reverse_descending_order[j][i] = i;
        }
        calculate_vector_descending_order(number_of_regimes + 1, m_log_Bq_reverse[j], m_log_Bq_reverse_descending_order[j]);
        // calculate log of sum Bq
        double largest_Bq_reverse = m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][0]];
        double temp_sum = 1;
        for (unsigned int index = 1; index < number_of_regimes + 1; index++){
            temp_sum += exp(m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][index]] - largest_Bq_reverse);
        }
        m_log_acceptance_prob -= largest_Bq_reverse + log(temp_sum);
    }
}

void rj::moving_full_changepoint_setup(const unsigned int & trace_index) {
    // choose which changepoint to move (we know that there is at least one that can be moved from our definitions of b_k, d_k, etc.)
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    m_log_proposal_ratio = 0;
    m_log_acceptance_prob = 0;
    // find the indices of the changepoints before and after m_tau[h] (so that we can sample the position of the new changepoint between them)
    unsigned long int tau_h_minus_1 = m_particle.get_changepoint(m_h - 1).get_position();
    unsigned long int tau_h_plus_1;
    if (m_h == m_dimension - 1){
        tau_h_plus_1 = m_end + 1;
    }
    else {
        tau_h_plus_1 = m_particle.get_changepoint(m_h + 1).get_position();
    }
    // sample the new changepoint position (and don't let m_tau[h] stay in the same position)
    m_new_changepoint_position = gsl_rng_uniform_int(r, tau_h_plus_1 - tau_h_minus_1 - 1) + tau_h_minus_1 + 1;
    if (tau_h_plus_1 - tau_h_minus_1 > 2) {
        while (m_particle.does_changepoint_exist_in_particle(m_new_changepoint_position, m_h, m_h)){
            m_new_changepoint_position = gsl_rng_uniform_int(r, tau_h_plus_1 - tau_h_minus_1 - 1) + tau_h_minus_1 + 1;
        }
    }
    m_log_k_prior_ratio = 0;
	m_adding_changepoint = changepoint(m_new_changepoint_position);
    
    // m_log_B_reverse[i] contains log Bayes factor if I_h_j_prime = i vs no changepoint here. m_log_q_reverse is the log I ratio for the same setting. m_log_Bq_reverse is their sum
    // m_log_B[i] contains log Bayes factor for adding a changepoint at m_new_changepoint_position with regime vs no changepoint there or at m_h.
    m_log_B = vector< vector< double > >(0);
    m_log_B_reverse = vector< vector< double > >(0);
    m_log_q = vector< vector< double > >(0);
    m_log_q_reverse = vector< vector< double > >(0);
    m_log_Bq = vector< vector< double > >(0);
    m_log_Bq_reverse = vector< vector< double > >(0);
    // m_log_Bq_reverse_descending_order[j][0] gives the index of the largest element in m_log_Bq_reverse, m_log_of_sum_Bq_reverse[j] gives the log of the sum of the Bq_reverse[j], same for m_log_Bq
    m_log_Bq_descending_order = vector< vector< unsigned int > >(0);
    m_log_Bq_reverse_descending_order = vector< vector< unsigned int > >(0);
    m_log_of_sum_Bq = vector< double >(m_number_of_processes, 0);
    m_log_of_sum_Bq_reverse = vector< double >(m_number_of_processes, 0);
    m_previous_log_likelihoods_without_right_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    m_log_likelihoods_with_right_sufficient_statistics = vector< vector< double > >(0);
    m_previous_log_likelihoods_with_right_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    m_previous_log_likelihoods_with_right_sufficient_statistics_reverse = vector< double >(m_number_of_processes, 0);
    m_actual_log_likelihoods_without_right_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    m_actual_log_likelihoods_without_right_sufficient_statistics_reverse = vector< double >(m_number_of_processes, 0);
    m_previous_log_likelihoods_without_middle_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    m_previous_log_likelihoods_with_middle_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    m_actual_log_likelihoods_with_middle_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    m_actual_log_likelihoods_without_middle_sufficient_statistics = vector< double >(m_number_of_processes, 0);
    // m_right_sufficient_statistics_reverse holds the sufficient statistics for the interval from m_h to m_h+1 for process j
    // the assignment here looks wrong, but it will have the cumulative data up to m_h subtracted later
    m_right_sufficient_statistics = vector< vector< double > >(m_number_of_processes);
    m_right_sufficient_statistics_reverse = vector< vector< double > >(m_number_of_processes);
    if (m_h == m_dimension - 1) {
        m_pm_ptr->get_cumulative_sufficient_data(tau_h_plus_1, m_right_sufficient_statistics);
        m_pm_ptr->get_cumulative_sufficient_data(tau_h_plus_1, m_right_sufficient_statistics_reverse);
    }
    else {
        m_pm_ptr->get_cumulative_sufficient_data(tau_h_plus_1, m_right_sufficient_statistics);
        m_pm_ptr->get_cumulative_sufficient_data(tau_h_plus_1, m_right_sufficient_statistics_reverse);
    }
	unsigned long int tau_h = m_particle.get_changepoint(m_h).get_position();
	m_tau_h_greater_than_tau_h_prime = tau_h > m_new_changepoint_position;
	m_middle_sufficient_statistics = vector< vector< double > >(m_number_of_processes);
	vector< vector< double > > sufficient_statistics_for_middle_sufficient_statistics(m_number_of_processes);
	if (m_tau_h_greater_than_tau_h_prime) {
		m_pm_ptr->get_cumulative_sufficient_data(tau_h, m_middle_sufficient_statistics);
		m_pm_ptr->get_cumulative_sufficient_data(m_new_changepoint_position, sufficient_statistics_for_middle_sufficient_statistics);
	}
	else {
		m_pm_ptr->get_cumulative_sufficient_data(m_new_changepoint_position, m_middle_sufficient_statistics);
		m_pm_ptr->get_cumulative_sufficient_data(tau_h, sufficient_statistics_for_middle_sufficient_statistics);
	}
    m_left_sufficient_statistics = vector< vector< double > >(m_number_of_processes);
    m_left_sufficient_statistics_reverse = vector< vector< double > >(m_number_of_processes);
    m_pm_ptr->get_cumulative_sufficient_data(m_new_changepoint_position, m_left_sufficient_statistics);
    m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h).get_position(), m_left_sufficient_statistics_reverse);
    vector< vector< double > > sufficient_statistics_up_to_left_cp_position(m_number_of_processes);
    m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h - 1).get_position(), sufficient_statistics_up_to_left_cp_position);
    m_removing_unobserved_regimes = vector< bool >(m_number_of_processes, false);
	for (unsigned int j = 0; j < m_number_of_processes; j++) {
		size_t number_of_regimes = m_particle.get_number_of_regimes(j);
		// get the actual regime for changepoint m_h
		unsigned int actual_regime = m_particle.get_previous_regime(m_h + 1, j);
        // if we are removing a whole regime from the process then we need to subtract 1 from the number of regimes.
        // if deleting tau_h means that regime r_j becomes unobserved then the number of regimes is reduced by 1.
        if (m_particle.removing_full_changepoint_leaves_highest_regime_unobserved(j, actual_regime)) {
            number_of_regimes--;
            m_removing_unobserved_regimes[j] = true;
        }
		// if we are removing a whole regime from the process then we need to subtract 1 from the number of regimes.
		m_log_B_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
		m_log_q_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
		m_log_B.push_back(vector< double >(number_of_regimes + 1, 0));
		m_log_q.push_back(vector< double >(number_of_regimes + 1, 0));
        m_log_likelihoods_with_right_sufficient_statistics.push_back(vector< double >(number_of_regimes + 1, 0));
		for (unsigned int index = 0; index < m_right_sufficient_statistics_reverse[j].size(); index++) {
			// this is correct, as m_left_sufficient_statistics_reverse currently holds the information up to m_h
			m_right_sufficient_statistics[j][index] -= m_left_sufficient_statistics[j][index];
			m_right_sufficient_statistics_reverse[j][index] -= m_left_sufficient_statistics_reverse[j][index];
			// we now make m_left_sufficient_statistics_reverse correct
			m_left_sufficient_statistics[j][index] -= sufficient_statistics_up_to_left_cp_position[j][index];
			m_left_sufficient_statistics_reverse[j][index] -= sufficient_statistics_up_to_left_cp_position[j][index];
			// also make m_middle_statistics correct
			m_middle_sufficient_statistics[j][index] -= sufficient_statistics_for_middle_sufficient_statistics[j][index];
		}
		// get the sufficient statistics for the regime that affects the changepoint prior to the new changepoint
		unsigned int previous_regime = m_particle.get_previous_regime(m_h, j);
		vector< double > previous_regime_sufficient_stats = m_particle.get_sufficient_statistics(j, previous_regime);
		// calculate the likelihood for the regime that affects the changepoint prior to the new changepoint
		double previous_regime_log_likelihood = m_particle.get_regime_log_likelihood(j, previous_regime);
		// calculate the likelihood for the previous and actual regime with the right sufficient statistics removed etc
		vector< double > previous_regime_sufficient_statistics_without_right_sufficient_statistics = previous_regime_sufficient_stats;
        vector< double > previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse = previous_regime_sufficient_stats;
		vector< double > previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse = previous_regime_sufficient_stats;
		vector< double > previous_regime_sufficient_statistics_with_middle_sufficient_statistics = previous_regime_sufficient_stats;
		vector< double > previous_regime_sufficient_statistics_without_middle_sufficient_statistics = previous_regime_sufficient_stats;
		for (unsigned int index = 0; index < previous_regime_sufficient_statistics_without_right_sufficient_statistics.size(); index++) {
            if (previous_regime == actual_regime) {
                previous_regime_sufficient_statistics_without_right_sufficient_statistics[index] -= m_right_sufficient_statistics[j][index];
                previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse[index] -= m_right_sufficient_statistics_reverse[j][index];
            }
            else {
                //previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse[index] -= m_right_sufficient_statistics_reverse[j][index]; delete?
                previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse[index] += m_right_sufficient_statistics_reverse[j][index];
                if (m_tau_h_greater_than_tau_h_prime) {
                    previous_regime_sufficient_statistics_without_middle_sufficient_statistics[index] -= m_middle_sufficient_statistics[j][index];
                }
                else {
                    previous_regime_sufficient_statistics_with_middle_sufficient_statistics[index] += m_middle_sufficient_statistics[j][index];
                }
            }
		}
        double previous_regime_log_likelihood_without_right_sufficient_statistics = 0, previous_regime_log_likelihood_without_right_sufficient_statistics_reverse = 0, previous_regime_log_likelihood_with_right_sufficient_statistics_reverse = 0, previous_regime_log_likelihood_with_middle_sufficient_statistics = 0, previous_regime_log_likelihood_without_middle_sufficient_statistics = 0;
        if (previous_regime == actual_regime) {
            previous_regime_log_likelihood_without_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_right_sufficient_statistics);
            previous_regime_log_likelihood_without_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse);
        }
        else {
            previous_regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse);
            if (m_tau_h_greater_than_tau_h_prime) {
                previous_regime_log_likelihood_without_middle_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_middle_sufficient_statistics);
            }
            else {
                previous_regime_log_likelihood_with_middle_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_with_middle_sufficient_statistics);
            }
        }
        
        vector< double > actual_regime_sufficient_statistics = vector< double >(0);
        vector< double > actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse = vector< double >(0);
        vector< double > actual_regime_sufficient_statistics_with_middle_sufficient_statistics = vector< double >(0);
        vector< double > actual_regime_sufficient_statistics_without_middle_sufficient_statistics = vector< double >(0);
        double actual_regime_log_likelihood = 0, actual_regime_log_likelihood_without_right_sufficient_statistics_reverse = 0, actual_regime_log_likelihood_with_middle_sufficient_statistics = 0, actual_regime_log_likelihood_without_middle_sufficient_statistics = 0;
        if (previous_regime != actual_regime) {
            actual_regime_sufficient_statistics = m_particle.get_sufficient_statistics(j, actual_regime);
            actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse = actual_regime_sufficient_statistics;
            if (m_tau_h_greater_than_tau_h_prime) {
                actual_regime_sufficient_statistics_with_middle_sufficient_statistics = actual_regime_sufficient_statistics;
            }
            else {
                actual_regime_sufficient_statistics_without_middle_sufficient_statistics = actual_regime_sufficient_statistics;
            }
            for (unsigned int index = 0; index < actual_regime_sufficient_statistics.size(); index++) {
                actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse[index] -= m_right_sufficient_statistics_reverse[j][index];
                if (m_tau_h_greater_than_tau_h_prime) {
                    actual_regime_sufficient_statistics_with_middle_sufficient_statistics[index] += m_middle_sufficient_statistics[j][index];
                }
                else {
                    actual_regime_sufficient_statistics_without_middle_sufficient_statistics[index] -= m_middle_sufficient_statistics[j][index];
                }
            }
            actual_regime_log_likelihood = m_pm_ptr->calculate_log_likelihood(j, actual_regime_sufficient_statistics);
            actual_regime_log_likelihood_without_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse);
            if (m_tau_h_greater_than_tau_h_prime) {
                actual_regime_log_likelihood_with_middle_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, actual_regime_sufficient_statistics_with_middle_sufficient_statistics);
            }
            else {
                actual_regime_log_likelihood_without_middle_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, actual_regime_sufficient_statistics_without_middle_sufficient_statistics);
            }
        }
        
		// for each regime calculate the reverse Bayes factor for assigning this regime to the new changepoint.
		for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
			if (previous_regime == actual_regime) {
                m_previous_log_likelihoods_without_right_sufficient_statistics[j] = previous_regime_log_likelihood_without_right_sufficient_statistics;
				if (regime == previous_regime) {
					m_log_B[j][regime] = 0;
					m_log_B_reverse[j][regime] = 0;
				}
				else {
					vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
					// calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
					vector< double > regime_sufficient_stats_with_right_sufficient_statistics = regime_sufficient_stats;
					vector< double > regime_sufficient_stats_with_right_sufficient_statistics_reverse = regime_sufficient_stats;
					for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
						regime_sufficient_stats_with_right_sufficient_statistics[index] += m_right_sufficient_statistics[j][index];
						regime_sufficient_stats_with_right_sufficient_statistics_reverse[index] += m_right_sufficient_statistics_reverse[j][index];
					}
                    double log_likelihoods_with_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics);
					m_log_B[j][regime] = log_likelihoods_with_right_sufficient_statistics + previous_regime_log_likelihood_without_right_sufficient_statistics - m_particle.get_regime_log_likelihood(j, regime) - previous_regime_log_likelihood;
					m_log_B_reverse[j][regime] = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse) + previous_regime_log_likelihood_without_right_sufficient_statistics_reverse - m_particle.get_regime_log_likelihood(j, regime) - previous_regime_log_likelihood;
                    m_log_likelihoods_with_right_sufficient_statistics[j][regime] = log_likelihoods_with_right_sufficient_statistics;
				}
			}
			else {
				if (regime == previous_regime) {
					m_log_B[j][regime] = 0;
					m_log_B_reverse[j][regime] = 0;
                    m_previous_log_likelihoods_with_right_sufficient_statistics_reverse[j] = previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                    m_actual_log_likelihoods_without_right_sufficient_statistics_reverse[j] = actual_regime_log_likelihood_without_right_sufficient_statistics_reverse;
				}
				else {
					if (regime == actual_regime) {
						if (m_tau_h_greater_than_tau_h_prime) {
							m_log_B[j][regime] = actual_regime_log_likelihood_with_middle_sufficient_statistics + previous_regime_log_likelihood_without_middle_sufficient_statistics - actual_regime_log_likelihood_without_right_sufficient_statistics_reverse - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
							m_log_B_reverse[j][regime] = actual_regime_log_likelihood + previous_regime_log_likelihood - actual_regime_log_likelihood_without_right_sufficient_statistics_reverse - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                            m_previous_log_likelihoods_without_middle_sufficient_statistics[j] = previous_regime_log_likelihood_without_middle_sufficient_statistics;
                            m_actual_log_likelihoods_with_middle_sufficient_statistics[j] = actual_regime_log_likelihood_with_middle_sufficient_statistics;
						}
						else {
							m_log_B[j][regime] = actual_regime_log_likelihood_without_middle_sufficient_statistics + previous_regime_log_likelihood_with_middle_sufficient_statistics - actual_regime_log_likelihood_without_right_sufficient_statistics_reverse - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
							m_log_B_reverse[j][regime] = actual_regime_log_likelihood + previous_regime_log_likelihood - actual_regime_log_likelihood_without_right_sufficient_statistics_reverse - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                            m_previous_log_likelihoods_with_middle_sufficient_statistics[j] = previous_regime_log_likelihood_with_middle_sufficient_statistics;
                            m_actual_log_likelihoods_without_middle_sufficient_statistics[j] = actual_regime_log_likelihood_without_middle_sufficient_statistics;
						}
					}
					else {
						vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
						vector< double > regime_sufficient_stats_with_right_sufficient_statistics = regime_sufficient_stats;
						vector< double > regime_sufficient_stats_with_right_sufficient_statistics_reverse = regime_sufficient_stats;
						for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
							regime_sufficient_stats_with_right_sufficient_statistics[index] += m_right_sufficient_statistics[j][index];
							regime_sufficient_stats_with_right_sufficient_statistics_reverse[index] += m_right_sufficient_statistics_reverse[j][index];
						}
                        double regime_log_likelihood_with_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics);
                        double regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse);
                        m_log_likelihoods_with_right_sufficient_statistics[j][regime] = regime_log_likelihood_with_right_sufficient_statistics;
                        m_actual_log_likelihoods_without_right_sufficient_statistics_reverse[j] = actual_regime_log_likelihood_without_right_sufficient_statistics_reverse;
						if (m_tau_h_greater_than_tau_h_prime) {
							m_log_B[j][regime] = regime_log_likelihood_with_right_sufficient_statistics + previous_regime_log_likelihood_without_middle_sufficient_statistics - m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats) - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
							m_log_B_reverse[j][regime] = regime_log_likelihood_with_right_sufficient_statistics_reverse + previous_regime_log_likelihood - m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats) - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                            m_previous_log_likelihoods_without_middle_sufficient_statistics[j] = previous_regime_log_likelihood_without_middle_sufficient_statistics;
						}
						else {
							m_log_B[j][regime] = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics) + previous_regime_log_likelihood_with_middle_sufficient_statistics - m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats) - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
							m_log_B_reverse[j][regime] = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse) + previous_regime_log_likelihood - m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats) - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                            m_previous_log_likelihoods_with_middle_sufficient_statistics[j] = previous_regime_log_likelihood_with_middle_sufficient_statistics;
						}
					}
				}
			}
		}
        // also propose making a new regime
        if (previous_regime == actual_regime) {
            vector< double > regime_sufficient_stats_with_right_sufficient_statistics = m_right_sufficient_statistics[j];
            vector< double > regime_sufficient_stats_with_right_sufficient_statistics_reverse = m_right_sufficient_statistics_reverse[j];
            double log_likelihoods_with_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics);
            m_log_B[j][number_of_regimes] = log_likelihoods_with_right_sufficient_statistics + previous_regime_log_likelihood_without_right_sufficient_statistics - previous_regime_log_likelihood;
            m_log_B_reverse[j][number_of_regimes] = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse) + previous_regime_log_likelihood_with_right_sufficient_statistics_reverse - previous_regime_log_likelihood;
            m_log_likelihoods_with_right_sufficient_statistics[j][number_of_regimes] = log_likelihoods_with_right_sufficient_statistics;
        }
        else {
            vector< double > regime_sufficient_stats_with_right_sufficient_statistics = m_right_sufficient_statistics[j];
            vector< double > regime_sufficient_stats_with_right_sufficient_statistics_reverse = m_right_sufficient_statistics_reverse[j];
            double regime_log_likelihood_with_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics);
            double regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse);
            m_log_likelihoods_with_right_sufficient_statistics[j][number_of_regimes] = regime_log_likelihood_with_right_sufficient_statistics;
            m_actual_log_likelihoods_without_right_sufficient_statistics_reverse[j] = actual_regime_log_likelihood_without_right_sufficient_statistics_reverse;
            if (m_tau_h_greater_than_tau_h_prime) {
                m_log_B[j][number_of_regimes] = regime_log_likelihood_with_right_sufficient_statistics + previous_regime_log_likelihood_without_middle_sufficient_statistics - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                m_log_B_reverse[j][number_of_regimes] = regime_log_likelihood_with_right_sufficient_statistics_reverse + previous_regime_log_likelihood - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                m_previous_log_likelihoods_without_middle_sufficient_statistics[j] = previous_regime_log_likelihood_without_middle_sufficient_statistics;
                if (actual_regime == number_of_regimes) { // because we have removed a regime and are now proposing adding it back in
                    m_actual_log_likelihoods_with_middle_sufficient_statistics[j] = actual_regime_log_likelihood_with_middle_sufficient_statistics;
                }
            }
            else {
                m_log_B[j][number_of_regimes] = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics) + previous_regime_log_likelihood_with_middle_sufficient_statistics - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                m_log_B_reverse[j][number_of_regimes] = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse) + previous_regime_log_likelihood - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                m_previous_log_likelihoods_with_middle_sufficient_statistics[j] = previous_regime_log_likelihood_with_middle_sufficient_statistics;
                if (actual_regime == number_of_regimes) { // because we have removed a regime and are now proposing adding it back in
                    m_actual_log_likelihoods_without_middle_sufficient_statistics[j] = actual_regime_log_likelihood_without_middle_sufficient_statistics;
                }
            }
        }

		if (m_h == m_dimension - 1 || m_particle.is_changepoint_index_separator_index(m_h + 1)) {
			for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
				m_log_q[j][regime] = m_log_q_reverse[j][regime] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, regime, false, actual_regime, m_removing_unobserved_regimes[j]);
			}
			m_log_q[j][number_of_regimes] = m_log_q_reverse[j][number_of_regimes] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, static_cast< unsigned int >(number_of_regimes), true, actual_regime, m_removing_unobserved_regimes[j]);
		}
		else {
			unsigned int subsequent_regime = m_particle.get_previous_regime(m_h + 2, j);
			for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
				m_log_q[j][regime] = m_log_q_reverse[j][regime] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, regime, false, actual_regime, m_removing_unobserved_regimes[j], subsequent_regime);
			}
			m_log_q[j][number_of_regimes] = m_log_q_reverse[j][number_of_regimes] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, static_cast<unsigned int>(number_of_regimes), true, actual_regime, m_removing_unobserved_regimes[j], subsequent_regime);
		}

		// set m_log_Bq = m_log_B + m_log_q
		m_log_Bq.push_back(vector< double >(number_of_regimes + 1, 0));
		m_log_Bq_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
		for (unsigned int regime = 0; regime < number_of_regimes + 1; regime++){
			m_log_Bq[j][regime] = m_log_B[j][regime] + m_log_q[j][regime];
			m_log_Bq_reverse[j][regime] = m_log_B_reverse[j][regime] + m_log_q_reverse[j][regime];
		}
        // check if regime number_of_regimes - 1 is unobserved. If removing tau_h reduces the number of regimes by 1, this will already be accounted for
        if (m_particle.is_regime_unobserved(j, number_of_regimes - 1)) {
            m_log_Bq[j][number_of_regimes - 1] = -1e300;
            m_log_Bq_reverse[j][number_of_regimes - 1] = -1e300;
        }
		// initialise m_log_Bq_reverse_descending_order and then calculate the ordering of the elements of m_log_Bq so that we can use fancy log and exp tricks to calculate the sum of the log Bq values.
		m_log_Bq_descending_order.push_back(vector< unsigned int >(number_of_regimes + 1, 0));
		m_log_Bq_reverse_descending_order.push_back(vector< unsigned int >(number_of_regimes + 1, 0));
		for (unsigned int i = 1; i < number_of_regimes + 1; i++) {
			m_log_Bq_descending_order[j][i] = i;
			m_log_Bq_reverse_descending_order[j][i] = i;
		}
		calculate_vector_descending_order(number_of_regimes + 1, m_log_Bq[j], m_log_Bq_descending_order[j]);
		calculate_vector_descending_order(number_of_regimes + 1, m_log_Bq_reverse[j], m_log_Bq_reverse_descending_order[j]);
		// calculate log of sum Bq
		double largest_Bq = m_log_Bq[j][m_log_Bq_descending_order[j][0]];
		double largest_Bq_reverse = m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][0]];
		double temp_sum = 1;
		double temp_sum_reverse = 1;
		for (unsigned int index = 1; index < number_of_regimes + 1; index++){
			temp_sum += exp(m_log_Bq[j][m_log_Bq_descending_order[j][index]] - largest_Bq);
			temp_sum_reverse += exp(m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][index]] - largest_Bq_reverse);
		}
        m_log_of_sum_Bq[j] = largest_Bq + log(temp_sum);
		m_log_acceptance_prob += largest_Bq - largest_Bq_reverse + log(temp_sum) - log(temp_sum_reverse);
	}
    // set the trace index
    m_particle.is_changepoint_index_separator_index(m_h); // need to run this because we ran is_changepoint_index_separator_index(m_h + 1); before and if m_h + 1 is a separator it will set the trace index to be one too high
}

void rj::resampling_full_changepoint_setup(const unsigned int & trace_index) {
	// choose which changepoint to resample
    // choose the changepoint index to resample
    int lower_index_bound = ((trace_index == 0) ? -1 : static_cast< int >(m_particle.get_separator_index(trace_index - 1)));
    unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
    m_h = static_cast< unsigned int >(gsl_rng_uniform_int(r, trace_dimension) + lower_index_bound + 1);
    m_log_proposal_ratio = 0;
	m_log_k_prior_ratio = 0;
	m_log_acceptance_prob = 0;

	// m_log_B_reverse[i] contains log Bayes factor if I_h_j_prime = i vs no changepoint here. m_log_q_reverse is the log I ratio for the same setting. m_log_Bq_reverse is their sum
	m_log_B_reverse = vector< vector< double > >(0);
	m_log_q_reverse = vector< vector< double > >(0);
	m_log_Bq_reverse = vector< vector< double > >(0);
    m_log_likelihoods_with_right_sufficient_statistics_reverse = vector< vector< double > >(0);
	// m_log_Bq_reverse_descending_order[j][0] gives the index of the largest element in m_log_Bq_reverse, m_log_of_sum_Bq_reverse[j] gives the log of the sum of the Bq_reverse[j]
	m_log_Bq_reverse_descending_order = vector< vector< unsigned int > >(0);
	m_log_of_sum_Bq_reverse = vector< double >(m_number_of_processes, 0);
    m_actual_log_likelihoods_without_right_sufficient_statistics_reverse = vector< double >(m_number_of_processes, 0);
	// m_right_sufficient_statistics_reverse holds the sufficient statistics for the interval from m_h to m_h+1 for process j
	// the assignment here looks wrong, but it will have the cumulative data up to m_h subtracted later
	m_right_sufficient_statistics_reverse = vector< vector< double > >(m_number_of_processes);
	// calculate how this move will affect the k prior
	if (m_h == m_dimension - 1) {
		m_pm_ptr->get_cumulative_sufficient_data(m_end + 1, m_right_sufficient_statistics_reverse);
	}
	else {
		m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h + 1).get_position(), m_right_sufficient_statistics_reverse);
	}
	m_left_sufficient_statistics_reverse = vector< vector< double > >(m_number_of_processes);
	m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h).get_position(), m_left_sufficient_statistics_reverse);
	vector< vector< double > > sufficient_statistics_up_to_left_cp_position(m_number_of_processes);
	m_pm_ptr->get_cumulative_sufficient_data(m_particle.get_changepoint(m_h - 1).get_position(), sufficient_statistics_up_to_left_cp_position);
    m_removing_unobserved_regimes = vector< bool >(m_number_of_processes, false);
	for (unsigned int j = 0; j < m_number_of_processes; j++) {
		size_t number_of_regimes = m_particle.get_number_of_regimes(j);
		// get the actual regime for changepoint m_h
		unsigned int actual_regime = m_particle.get_previous_regime(m_h + 1, j);
        // if we are removing a whole regime from the process then we need to subtract 1 from the number of regimes.
        // if deleting tau_h means that regime r_j becomes unobserved then the number of regimes is reduced by 1.
        if (m_particle.removing_full_changepoint_leaves_highest_regime_unobserved(j, actual_regime)) {
            number_of_regimes--;
            m_removing_unobserved_regimes[j] = true;
        }
		m_log_B_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
		m_log_q_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
        m_log_likelihoods_with_right_sufficient_statistics_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
		for (unsigned int index = 0; index < m_right_sufficient_statistics_reverse[j].size(); index++) {
			// this is correct, as m_left_sufficient_statistics_reverse currently holds the information up to m_h
			m_right_sufficient_statistics_reverse[j][index] -= m_left_sufficient_statistics_reverse[j][index];
			// we now make m_left_sufficient_statistics_reverse correct
			m_left_sufficient_statistics_reverse[j][index] -= sufficient_statistics_up_to_left_cp_position[j][index];
		}
		// get the sufficient statistics for the regime that affects the changepoint prior to the new changepoint
		unsigned int previous_regime = m_particle.get_previous_regime(m_h, j);
		vector< double > previous_regime_sufficient_stats = m_particle.get_sufficient_statistics(j, previous_regime);
		// calculate the likelihood for the regime that affects the changepoint prior to the new changepoint
		double previous_regime_log_likelihood = m_particle.get_regime_log_likelihood(j, previous_regime);
		// calculate the likelihood for the previous regime with the right sufficient statistics removed
		vector< double > previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse = previous_regime_sufficient_stats;
		for (unsigned int index = 0; index < previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse.size(); index++) {
			previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse[index] += m_right_sufficient_statistics_reverse[j][index];
		}
		double previous_regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse);
        
        vector< double > actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse = m_particle.get_sufficient_statistics(j, actual_regime);
        for (unsigned int index = 0; index < actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse.size(); index++) {
            actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse[index] -= m_right_sufficient_statistics_reverse[j][index];
        }
        m_actual_log_likelihoods_without_right_sufficient_statistics_reverse[j] = m_pm_ptr->calculate_log_likelihood(j, actual_regime_sufficient_statistics_without_right_sufficient_statistics_reverse);
		// for each regime calculate the reverse Bayes factor for assigning this regime to the new changepoint.
        for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
            if (previous_regime == actual_regime) {
                if (regime == previous_regime) {
                    m_log_B_reverse[j][regime] = 0;
                }
                else {
                    vector< double > previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse = previous_regime_sufficient_stats;
                    for (unsigned int index = 0; index < previous_regime_sufficient_statistics_with_right_sufficient_statistics_reverse.size(); index++) {
                        previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse[index] -= m_right_sufficient_statistics_reverse[j][index];
                    }
                    double previous_regime_log_likelihood_without_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse);
                    vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
                    // calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
                    vector< double > regime_sufficient_stats_with_right_sufficient_statistics_reverse = m_particle.get_sufficient_statistics(j, regime);
                    for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
                        regime_sufficient_stats_with_right_sufficient_statistics_reverse[index] += m_right_sufficient_statistics_reverse[j][index];
                    }
                    double regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse);
                    m_log_B_reverse[j][regime] = regime_log_likelihood_with_right_sufficient_statistics_reverse + previous_regime_log_likelihood_without_right_sufficient_statistics_reverse - m_particle.get_regime_log_likelihood(j, regime) - previous_regime_log_likelihood;
                    m_log_likelihoods_with_right_sufficient_statistics_reverse[j][regime] = regime_log_likelihood_with_right_sufficient_statistics_reverse;
                }
            }
            else {
                if (regime == previous_regime) {
                    m_log_B_reverse[j][regime] = 0;
                    m_log_likelihoods_with_right_sufficient_statistics_reverse[j][regime] = previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                }
                else {
                    if (regime == actual_regime) {
                        vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
                        // calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
                        vector< double > regime_sufficient_stats_without_right_sufficient_statistics = m_particle.get_sufficient_statistics(j, regime);
                        for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
                            regime_sufficient_stats_without_right_sufficient_statistics[index] -= m_right_sufficient_statistics_reverse[j][index];
                        }
                        m_log_B_reverse[j][regime] = m_particle.get_regime_log_likelihood(j, regime) + previous_regime_log_likelihood - m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_without_right_sufficient_statistics) - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
                    }
                    else {
                        vector< double > regime_sufficient_stats = m_particle.get_sufficient_statistics(j, regime);
                        // calculate the sufficient statistics for the regime if the right sufficient statistics are added to it
                        vector< double > regime_sufficient_stats_with_right_sufficient_statistics_reverse = m_particle.get_sufficient_statistics(j, regime);
                        for (unsigned int index = 0; index < regime_sufficient_stats.size(); index++) {
                            regime_sufficient_stats_with_right_sufficient_statistics_reverse[index] += m_right_sufficient_statistics_reverse[j][index];
                        }
                        double regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, regime_sufficient_stats_with_right_sufficient_statistics_reverse);
                        m_log_B_reverse[j][regime] = regime_log_likelihood_with_right_sufficient_statistics_reverse + previous_regime_log_likelihood - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse - m_particle.get_regime_log_likelihood(j, regime);
                        m_log_likelihoods_with_right_sufficient_statistics_reverse[j][regime] = regime_log_likelihood_with_right_sufficient_statistics_reverse;
                    }
                }
            }
        }
        // also propose adding a new regime
        if (previous_regime == actual_regime) {
            // calculate the likelihood for the previous regime with right interval removed
            vector< double > previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse = previous_regime_sufficient_stats;
            for (unsigned int index = 0; index < previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse.size(); index++) {
                previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse[index] -= m_right_sufficient_statistics_reverse[j][index];
            }
            double previous_regime_log_likelihood_without_right_sufficient_statistics = m_pm_ptr->calculate_log_likelihood(j, previous_regime_sufficient_statistics_without_right_sufficient_statistics_reverse);
            double regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, m_right_sufficient_statistics_reverse[j]);
            m_log_B_reverse[j][number_of_regimes] = previous_regime_log_likelihood_without_right_sufficient_statistics + regime_log_likelihood_with_right_sufficient_statistics_reverse - previous_regime_log_likelihood;
            m_log_likelihoods_with_right_sufficient_statistics_reverse[j][number_of_regimes] = regime_log_likelihood_with_right_sufficient_statistics_reverse;
        }
        else {
            double regime_log_likelihood_with_right_sufficient_statistics_reverse = m_pm_ptr->calculate_log_likelihood(j, m_right_sufficient_statistics_reverse[j]);
            m_log_B_reverse[j][number_of_regimes] = previous_regime_log_likelihood + regime_log_likelihood_with_right_sufficient_statistics_reverse - previous_regime_log_likelihood_with_right_sufficient_statistics_reverse;
            m_log_likelihoods_with_right_sufficient_statistics_reverse[j][number_of_regimes] = regime_log_likelihood_with_right_sufficient_statistics_reverse;
        }
        
        if (m_particle.is_changepoint_index_separator_index(m_h)) {
            if (m_particle.is_changepoint_index_separator_index(m_h + 1) || m_h == m_dimension - 1) {
                m_log_q_reverse[j][number_of_regimes] = m_particle.log_resampling_separator_changepoint_prior_ratio(j, true, true, m_removing_unobserved_regimes[j], actual_regime);
            }
            else {
                unsigned int following_regime = m_particle.get_previous_regime(m_h + 2, j);
                m_log_q_reverse[j][number_of_regimes] = m_particle.log_resampling_separator_changepoint_prior_ratio(j, false, true,  m_removing_unobserved_regimes[j], actual_regime, following_regime);
                for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                    m_log_q_reverse[j][regime] = m_particle.log_resampling_separator_changepoint_prior_ratio(j, false, false, m_removing_unobserved_regimes[j], actual_regime, following_regime, regime);
                }
            }
        }
        else {
            // for each regime calculate marked vector prior ratio
            if (m_h == m_dimension - 1 || m_particle.is_changepoint_index_separator_index(m_h + 1)) {
                for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                    m_log_q_reverse[j][regime] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, regime, false, actual_regime, m_removing_unobserved_regimes[j]);
                }
                m_log_q_reverse[j][number_of_regimes] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, static_cast< unsigned int >(number_of_regimes), true, actual_regime, m_removing_unobserved_regimes[j]);
            }
            else {
                unsigned int subsequent_regime = m_particle.get_previous_regime(m_h + 2, j);
                for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
                    m_log_q_reverse[j][regime] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, regime, false, actual_regime, m_removing_unobserved_regimes[j], subsequent_regime);
                }
                m_log_q_reverse[j][number_of_regimes] = m_particle.full_log_I_prior_remove_ratio(m_h, j, previous_regime, static_cast< unsigned int >(number_of_regimes), true, actual_regime, m_removing_unobserved_regimes[j], subsequent_regime);
            }
        }
        
        // set m_log_Bq = m_log_B + m_log_q
		m_log_Bq_reverse.push_back(vector< double >(number_of_regimes + 1, 0));
		for (unsigned int regime = 0; regime < number_of_regimes + 1; regime++){
			m_log_Bq_reverse[j][regime] = m_log_B_reverse[j][regime] + m_log_q_reverse[j][regime];
		}
        // check if regime number_of_regimes - 1 is unobserved. If removing tau_h reduces the number of regimes by 1, this will already be accounted for
        if (m_particle.is_regime_unobserved(j, number_of_regimes - 1)) {
            m_log_Bq_reverse[j][number_of_regimes - 1] = -1e300;
        }
		// initialise m_log_Bq_reverse_descending_order and then calculate the ordering of the elements of m_log_Bq so that we can use fancy log and exp tricks to calculate the sum of the log Bq values.
		m_log_Bq_reverse_descending_order.push_back(vector< unsigned int >(number_of_regimes + 1, 0));
		for (unsigned int i = 1; i < number_of_regimes + 1; i++) {
			m_log_Bq_reverse_descending_order[j][i] = i;
		}
		calculate_vector_descending_order(number_of_regimes + 1, m_log_Bq_reverse[j], m_log_Bq_reverse_descending_order[j]);
		// calculate log of sum Bq
		double largest_Bq_reverse = m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][0]];
		double temp_sum = 1;
		for (unsigned int index = 1; index < number_of_regimes + 1; index++){
			temp_sum += exp(m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][index]] - largest_Bq_reverse);
		}
        m_log_of_sum_Bq_reverse[j] = largest_Bq_reverse + log(temp_sum);
	}
}

void rj::altering_unobserved_regimes_setup() {
    m_log_acceptance_prob = 0;
    m_log_k_prior_ratio = 0;
    m_log_regimes_prior_ratio = 0;
    m_log_full_I_prior_ratio = 0;
    m_altering_unobserved_regimes = vector< int >(m_number_of_processes, 0);
    double alpha, beta, u, log_acceptance_prob, log_full_I_prior_ratio, log_regimes_prior_ratio;
    for (unsigned int process = 0; process < m_number_of_processes; process++) {
        if (m_particle.get_number_of_unobserved_regimes()[process] > 0) {
            alpha = 1.0 / 2.0, beta = 1.0;
        }
        else {
            alpha = 1.0, beta = 1.0;
        }
        u = gsl_ran_flat(r, 0, 1);
        double r_j = static_cast< double >(m_particle.get_number_of_regimes(process));
        double r_tilde_j = static_cast< double >(m_particle.get_number_of_unobserved_regimes()[process]);
        if (u < alpha) {
            log_acceptance_prob = log(0.5) - log(alpha) + log(r_j) - log(r_tilde_j + 1.0);
            log_full_I_prior_ratio = m_particle.calculate_and_get_add_unobserved_regimes_full_I_prior_ratio(process);
            log_regimes_prior_ratio = log(1.0 - m_particle.get_rho());
            log_acceptance_prob += log_full_I_prior_ratio + log_regimes_prior_ratio;
            if (log_acceptance_prob > 0 || (log(gsl_ran_flat(r, 0, 1)) < log_acceptance_prob)) {
                m_altering_unobserved_regimes[process] = 1;
                m_log_full_I_prior_ratio += log_full_I_prior_ratio;
                m_log_regimes_prior_ratio += log_regimes_prior_ratio;
            }
        } else if (u < beta) {
            log_acceptance_prob = log(0.5) - log(beta - alpha) + (m_particle.get_number_of_unobserved_regimes()[process] == 1 ? log(2.0) : 0.0) + log(r_tilde_j) - log(r_j - 1.0);
            log_full_I_prior_ratio = m_particle.calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio(process);
            log_regimes_prior_ratio = -log(1 - m_particle.get_rho());
            log_acceptance_prob += log_full_I_prior_ratio + log_regimes_prior_ratio;
            if (log_acceptance_prob > 0 || (log(gsl_ran_flat(r, 0, 1)) < log_acceptance_prob)) {
                m_altering_unobserved_regimes[process] = -1;
                m_log_full_I_prior_ratio += log_full_I_prior_ratio;
                m_log_regimes_prior_ratio += log_regimes_prior_ratio;
            }
        }
    }
}

void rj::full_acceptance_procedure(const double & u1){
    m_particle.increase_log_k_prior(m_log_k_prior_ratio);
    if (u1 < m_b_k) {
        // choose the new regimes for this new changepoint
        vector< unsigned int > new_regimes(m_number_of_processes, 0);
        double temp_sum, regime_chooser;
        unsigned int index;
        bool chosen;
        for (unsigned int j = 0; j < m_number_of_processes; j++) {
            temp_sum = 0;
            regime_chooser = log(gsl_ran_flat(r, 0, 1)) + m_log_of_sum_Bq[j];
            index = 0;
            chosen = false;
            do {
                temp_sum += exp(m_log_Bq[j][m_log_Bq_descending_order[j][index]] - m_log_Bq[j][m_log_Bq_descending_order[j][0]]);
                chosen = regime_chooser <= m_log_Bq[j][m_log_Bq_descending_order[j][0]] + log(temp_sum);
                index++;
            }
            while(!chosen);
            new_regimes[j] = m_log_Bq_descending_order[j][index - 1];
            m_particle.increase_log_likelihood(m_log_B[j][new_regimes[j]]);
            bool adding_new_regime = new_regimes[j] == m_log_q[j].size() - 1;
            m_particle.increase_log_full_I_prior(m_log_q[j][new_regimes[j]], adding_new_regime, j); // if a new regime is added, the (1-rho) in m_log_q[j][new_regimes[j]] will be added to m_log_regime_prior and taken away from m_log_full_I_prior
        }
        vector< double > right_number_of_observations = vector< double >(m_number_of_processes);
        m_pm_ptr->get_number_of_observations(m_right_sufficient_statistics, right_number_of_observations);
        m_particle.add_full_changepoint(m_particle.get_add_changepoint_index(), m_adding_changepoint, new_regimes, m_right_sufficient_statistics, m_log_likelihoods_with_right_sufficient_statistics, m_previous_log_likelihoods_without_right_sufficient_statistics, right_number_of_observations);
    }
    else if (u1 < m_d_k) {
        for (unsigned int j = 0; j < m_number_of_processes; j++) {
            unsigned int actual_regime = m_particle.get_previous_regime(m_h + 1, j);
            m_particle.increase_log_likelihood(-m_log_B_reverse[j][actual_regime]);
            m_particle.decrease_log_full_I_prior(m_log_q_reverse[j][actual_regime], m_removing_unobserved_regimes[j], false, j);
        }
        vector< double > right_number_of_observations_reverse = vector< double >(m_number_of_processes);
        m_pm_ptr->get_number_of_observations(m_right_sufficient_statistics_reverse, right_number_of_observations_reverse);
        m_particle.remove_full_changepoint(m_h, m_right_sufficient_statistics_reverse, m_actual_log_likelihoods_without_right_sufficient_statistics_reverse, m_previous_log_likelihoods_with_right_sufficient_statistics_reverse, right_number_of_observations_reverse, m_removing_unobserved_regimes);
    }
    else if (u1 < m_m_k) {
        // choose the new regimes for this new changepoint
        vector< unsigned int > new_regimes(m_number_of_processes, 0);
        double temp_sum, regime_chooser;
        unsigned int index;
        bool chosen;
        for (unsigned int j = 0; j < m_number_of_processes; j++) {
            temp_sum = 0;
            regime_chooser = log(gsl_ran_flat(r, 0, 1)) + m_log_of_sum_Bq[j];
            index = 0;
            chosen = false;
            do {
                temp_sum += exp(m_log_Bq[j][m_log_Bq_descending_order[j][index]] - m_log_Bq[j][m_log_Bq_descending_order[j][0]]);
                chosen = regime_chooser <= m_log_Bq[j][m_log_Bq_descending_order[j][0]] + log(temp_sum);
                index++;
            }
            while(!chosen);
            new_regimes[j] = m_log_Bq_descending_order[j][index - 1];
            unsigned int actual_regime = m_particle.get_previous_regime(m_h + 1, j);
            m_particle.increase_log_likelihood(m_log_B[j][new_regimes[j]] - m_log_B_reverse[j][actual_regime]);
            bool adding_new_regime = new_regimes[j] == m_log_q[j].size() - 1;
            m_particle.decrease_log_full_I_prior(-m_log_q[j][new_regimes[j]] + m_log_q_reverse[j][actual_regime], m_removing_unobserved_regimes[j], adding_new_regime, j); // if a new regime is added, the (1-rho) in m_log_q[j][new_regimes[j]] will be added to m_log_regime_prior and taken away from m_log_full_I_prior
        }
        vector< double > right_number_of_observations = vector< double >(m_number_of_processes);
        m_pm_ptr->get_number_of_observations(m_right_sufficient_statistics, right_number_of_observations);
        vector< double > right_number_of_observations_reverse = vector< double >(m_number_of_processes);
        m_pm_ptr->get_number_of_observations(m_right_sufficient_statistics_reverse, right_number_of_observations_reverse);
        vector< double > middle_number_of_observations = vector< double >(m_number_of_processes);
        m_pm_ptr->get_number_of_observations(m_middle_sufficient_statistics, middle_number_of_observations);
		m_particle.move_full_changepoint(m_h, m_adding_changepoint, new_regimes, m_tau_h_greater_than_tau_h_prime, m_right_sufficient_statistics, right_number_of_observations, m_right_sufficient_statistics_reverse, right_number_of_observations_reverse, m_middle_sufficient_statistics, middle_number_of_observations,m_previous_log_likelihoods_without_right_sufficient_statistics, m_previous_log_likelihoods_with_right_sufficient_statistics_reverse, m_previous_log_likelihoods_with_middle_sufficient_statistics, m_previous_log_likelihoods_without_middle_sufficient_statistics, m_log_likelihoods_with_right_sufficient_statistics, m_actual_log_likelihoods_with_middle_sufficient_statistics, m_actual_log_likelihoods_without_middle_sufficient_statistics, m_actual_log_likelihoods_without_right_sufficient_statistics_reverse, m_removing_unobserved_regimes);
    }
    else if (u1 < m_r_k){
        // choose the new regimes for this new changepoint
        vector< unsigned int > new_regimes(m_number_of_processes, 0);
        double temp_sum, regime_chooser;
        unsigned int index;
        bool chosen;
        for (unsigned int j = 0; j < m_number_of_processes; j++) {
            temp_sum = 0;
            regime_chooser = log(gsl_ran_flat(r, 0, 1)) + m_log_of_sum_Bq_reverse[j];
            index = 0;
            chosen = false;
            do {
                temp_sum += exp(m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][index]] - m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][0]]);
                chosen = regime_chooser <= m_log_Bq_reverse[j][m_log_Bq_reverse_descending_order[j][0]] + log(temp_sum);
                index++;
            }
            while(!chosen);
            new_regimes[j] = m_log_Bq_reverse_descending_order[j][index - 1];
            unsigned int actual_regime = m_particle.get_previous_regime(m_h + 1, j);
            m_particle.increase_log_likelihood(m_log_B_reverse[j][new_regimes[j]] - m_log_B_reverse[j][actual_regime]);
            bool adding_new_regime = new_regimes[j] == m_log_q_reverse[j].size() - 1;
            m_particle.decrease_log_full_I_prior(-m_log_q_reverse[j][new_regimes[j]] + m_log_q_reverse[j][actual_regime], m_removing_unobserved_regimes[j], adding_new_regime, j); // if a new regime is added, the (1-rho) in m_log_q[j][new_regimes[j]] will be added to m_log_regime_prior and taken away from m_log_full_I_prior
        }
        vector< double > right_number_of_observations_reverse = vector< double >(m_number_of_processes);
        m_pm_ptr->get_number_of_observations(m_right_sufficient_statistics_reverse, right_number_of_observations_reverse);
        m_particle.resample_full_changepoint(m_h, new_regimes, m_right_sufficient_statistics_reverse, right_number_of_observations_reverse, m_log_likelihoods_with_right_sufficient_statistics_reverse, m_actual_log_likelihoods_without_right_sufficient_statistics_reverse, m_removing_unobserved_regimes);
    }
    else if (u1 < m_au_k) {
        m_particle.increase_log_full_I_prior_unobserved(m_log_full_I_prior_ratio, m_altering_unobserved_regimes);
        m_particle.increase_log_regimes_prior(m_log_regimes_prior_ratio);
        m_particle.alter_unobserved_regimes(m_altering_unobserved_regimes, m_number_of_processes);
    }
}

void rj::full_recording_procedure() {
	m_recorded_full_dimensions.push_back(m_particle.get_dimension());
    //calculate the effective dimension (i.e. don't count changepoints where none of the processes change regime at that changepoint.
    m_recorded_full_effective_dimensions.push_back(m_particle.calculate_and_get_full_effective_dimension());
    
	vector< unsigned long int > changepoint_hist = m_particle.calculate_and_get_full_changepoint_histogram(m_number_of_changepoint_bins, m_number_of_processes);
	size_t size_of_recorded_changepoints = m_recorded_full_changepoints.size();
	if (size_of_recorded_changepoints > 0){ //have we started recording changepoint histograms, or is this the first time?
		for (unsigned long int i = 0; i < size_of_recorded_changepoints; i++){
			m_recorded_full_changepoints[i] += changepoint_hist[i];
		}
	}
	else {
		m_recorded_full_changepoints = changepoint_hist;
	}
    
    vector< size_t > number_of_regimes(m_number_of_processes, 0);
    for (unsigned int proc = 0; proc < m_number_of_processes; proc++) {
        number_of_regimes[proc] = m_particle.get_number_of_regimes(proc);
    }
    m_recorded_number_of_regimes.push_back(number_of_regimes);
    vector< unsigned int > number_of_observed_regimes(m_number_of_processes, 0);
    for (unsigned int proc = 0; proc < m_number_of_processes; proc++) {
        number_of_observed_regimes[proc] = m_particle.get_number_of_observed_regimes(proc);
    }
    m_recorded_number_of_observed_regimes.push_back(number_of_observed_regimes);
    
    if (m_number_of_traces > 0) {
        if (m_recorded_similarity_matrix.size() > 0) {
            m_particle.calculate_and_add_similarity_matrices(m_recorded_similarity_matrix, m_number_of_processes);
            m_particle.calculate_and_add_min_proportion_similarity_matrices(m_recorded_min_proportion_similarity_matrix, m_number_of_processes, m_observations_in_each_trace);
			m_recorded_similarity_matrices.push_back(vector< vector< double > >(m_number_of_traces, vector< double >(m_number_of_traces, 0.0)));
			m_recorded_min_proportion_similarity_matrices.push_back(vector< vector< double > >(m_number_of_traces, vector< double >(m_number_of_traces, 0.0)));
			m_particle.calculate_and_add_similarity_matrices(m_recorded_similarity_matrices.back(), m_number_of_processes);
			m_particle.calculate_and_add_min_proportion_similarity_matrices(m_recorded_min_proportion_similarity_matrices.back(), m_number_of_processes, m_observations_in_each_trace);
        }
        else {
            m_recorded_similarity_matrix = vector< vector< double > >(m_number_of_traces, vector< double >(m_number_of_traces, 0.0));
            m_recorded_min_proportion_similarity_matrix = vector< vector< double > >(m_number_of_traces, vector< double >(m_number_of_traces, 0.0));
			m_recorded_similarity_matrices = vector< vector< vector< double > > >(1, vector< vector< double > >(m_number_of_traces, vector< double >(m_number_of_traces, 0.0)));
			m_recorded_min_proportion_similarity_matrices = vector< vector< vector< double > > >(1, vector< vector< double > >(m_number_of_traces, vector< double >(m_number_of_traces, 0.0)));
            m_particle.calculate_and_add_similarity_matrices(m_recorded_similarity_matrix, m_number_of_processes);
            m_particle.calculate_and_add_min_proportion_similarity_matrices(m_recorded_min_proportion_similarity_matrix, m_number_of_processes, m_observations_in_each_trace);
			m_particle.calculate_and_add_similarity_matrices(m_recorded_similarity_matrices.back(), m_number_of_processes);
			m_particle.calculate_and_add_min_proportion_similarity_matrices(m_recorded_min_proportion_similarity_matrices.back(), m_number_of_processes, m_observations_in_each_trace);
        }
    }
    
    if (m_recording_association_matrix) {
        for (unsigned int process = 0; process < m_number_of_processes; process++) {
            m_particle.add_to_association_matrix(m_association_matrices[process], process);
        }
    }
    
	double log_posterior = m_particle.get_full_log_posterior();
	m_recorded_full_log_posteriors.push_back(log_posterior);
}

void rj::update_full_MAP() {
    double log_posterior = m_particle.get_full_log_posterior();
    if (m_full_MAP_log_posterior < log_posterior){
        m_full_MAP_particle = m_particle;
        m_full_MAP_dimension = m_full_MAP_particle.get_dimension();
        m_full_MAP_log_posterior = log_posterior;
    }
}

void rj::run_full_simulation(){
    cout << "starting full changepoint stage" << endl;
    m_particle.print_likelihood();
    cout << "the posterior is " << m_particle.get_full_log_posterior() << endl << endl;
    // burnin
    for (unsigned long int iteration = 0; iteration < m_full_burnin; iteration++) {
        m_dimension = m_particle.get_dimension();
        unsigned int trace_index = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_number_of_traces));
        unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
        if (trace_dimension > 0) {
            m_b_k = 1.0 / 4.0, m_d_k = 2.0 / 4.0, m_m_k = 2.75 / 4.0, m_r_k = 3.5 / 4.0, m_au_k = 4.0 / 4.0;
        } else {
            m_b_k = 1.0;
        }
        
        double u1 = gsl_ran_flat(r, 0, 1);
        if (u1 < m_b_k) { // birth
            adding_full_changepoint_setup(trace_index);
        } else if (u1 < m_d_k) { // death
            removing_full_changepoint_setup(trace_index);
        } else if (u1 < m_m_k) { // move
            moving_full_changepoint_setup(trace_index);
        } else if (u1 < m_r_k) { // resample marked vector
            resampling_full_changepoint_setup(trace_index);
        } else if (u1 < m_au_k) { // add unobserved regimes
            altering_unobserved_regimes_setup();
        }

        if (m_log_acceptance_prob >= 0 || (log(gsl_ran_flat(r, 0, 1)) < m_log_acceptance_prob)){
            full_acceptance_procedure(u1);
        }
        
        //m_particle.check_observations_in_traces(m_end + 1);
        //check_total_full_log_likelihood(m_particle);
        //m_particle.check_full_log_posterior();
        /*if (iteration % 1000 == 0 && (m_particle.calculate_and_get_full_log_posterior(m_number_of_processes) - m_particle.get_full_log_posterior() > 0.000001 || m_particle.calculate_and_get_full_log_posterior(m_number_of_processes) - m_particle.get_full_log_posterior() < -0.000001)) {
            cout << iteration << '\t' << m_particle.calculate_and_get_full_log_posterior(m_number_of_processes) - m_particle.get_full_log_posterior() << endl;
        }*/
    }
    
    for (unsigned long int iteration = 0; iteration < m_full_iterations; iteration++) {
        m_dimension = m_particle.get_dimension();
        unsigned int trace_index = static_cast< unsigned int >(gsl_rng_uniform_int(r, m_number_of_traces));
        unsigned int trace_dimension = m_particle.calculate_trace_dimension(trace_index);
        if (trace_dimension > 0) {
            m_b_k = 1.0 / 4.0, m_d_k = 2.0 / 4.0, m_m_k = 2.75 / 4.0, m_r_k = 3.5 / 4.0, m_au_k = 4.0 / 4.0;
        } else {
            m_b_k = 1.0;
        }
        
        double u1 = gsl_ran_flat(r, 0, 1);
        if (u1 < m_b_k) { // birth
            adding_full_changepoint_setup(trace_index);
            m_recorded_full_birth_proposals++;
        } else if (u1 < m_d_k) { // death
            removing_full_changepoint_setup(trace_index);
            m_recorded_full_death_proposals++;
        } else if (u1 < m_m_k) { // move
            moving_full_changepoint_setup(trace_index);
            m_recorded_full_move_proposals++;
        } else if (u1 < m_r_k) { // resample marked vector
            resampling_full_changepoint_setup(trace_index);
            m_recorded_full_resample_proposals++;
        } else if (u1 < m_au_k) { // add unobserved regimes
            altering_unobserved_regimes_setup();
            m_recorded_full_unobserveds_proposals++;
        }

        if (m_log_acceptance_prob > 0 || (log(gsl_ran_flat(r, 0, 1)) < m_log_acceptance_prob)){
            full_acceptance_procedure(u1);
            // check if this alters the best-seen set of change points and regimes
            update_full_MAP();
            if (u1 < m_b_k) { // birth
                m_recorded_full_birth_acceptances++;
            } else if (u1 < m_d_k) { // death
                m_recorded_full_death_acceptances++;
            } else if (u1 < m_m_k) { // move
                m_recorded_full_move_acceptances++;
            } else if (u1 < m_r_k) { // resample marked vector
                m_recorded_full_resample_acceptances++;
            } else if (u1 < m_au_k) { // add unobserved regimes
                m_recorded_full_unobserveds_acceptances++;
            }
        }
        //check_total_full_log_likelihood(m_particle);
        //m_particle.check_full_log_posterior();
        if (m_recording_full_samples && (iteration % m_full_thinning == 0)) { //store sample
            full_recording_procedure();
            /* CODE IF YOU WANT TO REGULARLY WRITE OUTPUT TO FILE
             if (CONDITION ON ITERATION) {
             string MAP_cps_Filename = m_data_file + "_full_MAP_CPs.txt";
             rjobject.write_full_MAP_changepoints_to_file(MAP_cps_Filename);
             string dimension_distribution_Filename = m_data_file + "_full_dimension_distribution.txt";
             rjobject.write_full_dimension_distribution_to_file(dimension_distribution_Filename);
             string changepoints_distribution_Filename = m_data_file + "_full_changepoints_distribution.txt";
             rjobject.write_full_changepoints_distribution_to_file(changepoints_distribution_Filename, iteration);
             string log_posterior_trace_Filename = m_data_file + "_full_log_posterior_trace.txt";
             rjobject.write_full_log_posterior_trace_to_file(log_posterior_trace_Filename);
             }*/
        }
        
        
        /*if (m_particle.calculate_and_get_full_log_posterior(m_number_of_processes) - m_particle.get_full_log_posterior() > 0.000001 || m_particle.calculate_and_get_full_log_posterior(m_number_of_processes) - m_particle.get_full_log_posterior() < -0.000001) {
            cout << iteration << '\t' << m_particle.calculate_and_get_full_log_posterior(m_number_of_processes) - m_particle.get_full_log_posterior() << endl;
        }*/
    }
    cout << "ending full changepoint stage" << endl;
    m_particle.print_likelihood();
    cout << "the posterior is " << m_particle.get_full_log_posterior() << endl;
}

/*void rj::add_unobserved_regimes_setup() {
 m_adding_unobserved_regimes = vector< bool >(m_number_of_processes, false);
 for (unsigned int process = 0; process < m_number_of_processes; process++) {
 m_adding_unobserved_regimes[process] = gsl_ran_flat(r, 0, 1) > 0.5;
 }
 m_log_k_prior_ratio = 0;
 m_log_full_I_prior_ratio = m_particle.calculate_and_get_add_unobserved_regimes_full_I_prior_ratio(m_adding_unobserved_regimes, m_number_of_processes);
 m_log_regimes_prior_ratio = m_particle.calculate_and_get_add_unobserved_regimes_regimes_prior_ratio(m_adding_unobserved_regimes, m_number_of_processes);
 m_log_acceptance_prob = m_log_full_I_prior_ratio + m_log_regimes_prior_ratio + m_particle.calculate_and_get_add_unobserved_regimes_proposal_ratio(m_adding_unobserved_regimes, m_number_of_processes);
 }
 
 void rj::remove_unobserved_regimes_setup() {
 vector< bool > any_unobserved_regimes = m_particle.any_unobserved_regimes();
 m_removing_unobserved_regimes = vector< bool >(m_number_of_processes, false);
 for (unsigned int process = 0; process < m_number_of_processes; process++) {
 if (any_unobserved_regimes[process]) {
 m_removing_unobserved_regimes[process] = gsl_ran_flat(r, 0, 1) > 0.5;
 }
 }
 m_log_k_prior_ratio = 0;
 m_log_full_I_prior_ratio = m_particle.calculate_and_get_remove_unobserved_regimes_full_I_prior_ratio(m_removing_unobserved_regimes, m_number_of_processes);
 m_log_regimes_prior_ratio = m_particle.calculate_and_get_remove_unobserved_regimes_regimes_prior_ratio(m_removing_unobserved_regimes, m_number_of_processes);
 m_log_acceptance_prob = m_log_full_I_prior_ratio + m_log_regimes_prior_ratio + m_particle.calculate_and_get_remove_unobserved_regimes_proposal_ratio(m_removing_unobserved_regimes, m_number_of_processes);
 }*/
/*else if (u1 < m_au_k) {
 m_particle.increase_log_full_I_prior_unobserved(m_log_full_I_prior_ratio, m_adding_unobserved_regimes, true);
 m_particle.increase_log_regimes_prior(m_log_regimes_prior_ratio);
 m_particle.add_unobserved_regimes(m_adding_unobserved_regimes, m_number_of_processes);
 }
 else if (u1 < m_ru_k) {
 m_particle.increase_log_full_I_prior_unobserved(m_log_full_I_prior_ratio, m_removing_unobserved_regimes, false);
 m_particle.increase_log_regimes_prior(m_log_regimes_prior_ratio);
 m_particle.remove_unobserved_regimes(m_removing_unobserved_regimes, m_number_of_processes);
 }*/
/*else if (u1 < m_au_k) {
 for (unsigned int process = 0; process < m_number_of_processes; process++) {
 if (m_adding_unobserved_regimes[process]) {
 full_log_acceptance_probability += log(1 - m_particle.get_rho());
 size_t number_of_regimes = m_particle.get_number_of_regimes(process);
 double r_j = static_cast< double >(number_of_regimes), delta = m_particle.get_dirichlet_alpha();
 for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
 double k_j_s = static_cast< double >(m_particle.get_number_of_right_transitions(process, regime));
 full_log_acceptance_probability += gsl_sf_lngamma(r_j * delta + delta) + gsl_sf_lngamma(r_j * delta + k_j_s) - gsl_sf_lngamma(r_j * delta) - gsl_sf_lngamma(r_j * delta + delta + k_j_s);
 }
 double r_tilde_j = static_cast< double >(m_particle.get_number_of_unobserved_regimes()[process]);
 full_log_acceptance_probability += log(r_j) - log(r_tilde_j + 1);
 }
 }
 }
 else if (u1 < m_ru_k) {
 for (unsigned int process = 0; process < m_number_of_processes; process++) {
 if (m_removing_unobserved_regimes[process]) {
 full_log_acceptance_probability -= log(1 - m_particle.get_rho());
 size_t number_of_regimes = m_particle.get_number_of_regimes(process);
 double r_j = static_cast< double >(number_of_regimes), delta = m_particle.get_dirichlet_alpha();
 for (unsigned int regime = 0; regime < number_of_regimes; regime++) {
 double k_j_s = static_cast< double >(m_particle.get_number_of_right_transitions(process, regime));
 full_log_acceptance_probability += gsl_sf_lngamma(r_j * delta - delta) + gsl_sf_lngamma(r_j * delta + k_j_s) - gsl_sf_lngamma(r_j * delta) - gsl_sf_lngamma(r_j * delta - delta + k_j_s);
 }
 double r_tilde_j = static_cast< double >(m_particle.get_number_of_unobserved_regimes()[process]);
 full_log_acceptance_probability += log(r_tilde_j) - log(r_j - 1);
 }
 */

#endif
