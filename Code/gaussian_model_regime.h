#ifndef GAUSSIAN_MODEL_REGIME_h
#define GAUSSIAN_MODEL_REGIME_h

#include "basic_probability_model.h"
#include "data_input.h"

class gaussian_model: public probability_model{
    
public:
    gaussian_model(const vector< string > & filenames, const unsigned long int & end, const unsigned long int & diff, const double & mu_0, const double & kappa_0, const double & alpha_0, const double & beta_0);
    ~gaussian_model();
    double calculate_log_likelihood(const unsigned long int & position1, const unsigned long int & position2);
    double calculate_log_likelihood(const vector< double > & sufficient_statistics);
    void print_data();
    void get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations); // calculate the number of observations given a set of sufficient statistics
    
protected:
    double m_mu_0;
    double m_kappa_0;
    double m_alpha_0;
    double m_beta_0;
};

gaussian_model::gaussian_model(const vector< string > & filenames, const unsigned long int & end, const unsigned long int & diff, const double & mu_0, const double & kappa_0, const double & alpha_0, const double & beta_0):probability_model(end, diff), m_mu_0(mu_0), m_kappa_0(kappa_0), m_alpha_0(alpha_0), m_beta_0(beta_0) {
    m_data = new data_input(filenames[0], filenames[1], end, "gaussian", 3, diff);
}

gaussian_model::~gaussian_model() {
    m_data->~data_input();
}

double gaussian_model::calculate_log_likelihood(const unsigned long int & position1, const unsigned long int & position2) {
    vector< double > sufficient_stats1;
    m_data->get_cumulative_data_row(position1, sufficient_stats1);
    vector< double > sufficient_stats2;
    m_data->get_cumulative_data_row(position2, sufficient_stats2);
    for (size_t i = 0; i < sufficient_stats1.size(); i++) {
        sufficient_stats2[i] -= sufficient_stats1[i];
    }
    return calculate_log_likelihood(sufficient_stats2);
}

double gaussian_model::calculate_log_likelihood(const vector< double > & sufficient_stats) { // using notation from http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
    double n = sufficient_stats[2];
    if (n == 0) { // check if the interval is empty
        return 0.0;
    }
    double x_bar = sufficient_stats[0] / n;
    double sum_of_x_i_minus_x_bar_squared = sufficient_stats[1] - n * x_bar * x_bar;
    double pi = 3.14159265358979323846;
    double kappa_n = m_kappa_0 + n;
    double alpha_n = m_alpha_0 + n / 2.0;
    double beta_n = m_beta_0 + sum_of_x_i_minus_x_bar_squared / 2.0 + (m_kappa_0 * n * (x_bar - m_mu_0) * (x_bar - m_mu_0)) / (2.0 * (m_kappa_0 + n));
    return gsl_sf_lngamma(alpha_n) - gsl_sf_lngamma(m_alpha_0) + m_alpha_0 * log(m_beta_0) - alpha_n * log(beta_n) + 0.5 * log(m_kappa_0 / kappa_n) - (n / 2.0) * log(2 * pi);
}

void gaussian_model::print_data() {
    m_data->print_data();
}

void gaussian_model::get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations) {
    number_of_observations = sufficient_stats[2];
}

#endif
