#ifndef MARKOV_CHAIN_MODEL_H
#define MARKOV_CHAIN_MODEL_H

#include <math.h>
//#include "probability_model.h"
#include <stdlib.h>
using namespace std;

class markov_chain : public probability_model{

public:
    markov_chain( bool independent, vector< string > filenames, unsigned long int number_of_states, unsigned long int end, unsigned long int diff, double markov_alpha, vector< unsigned long int > separators = vector< unsigned long int >(0));
    ~markov_chain();
    double calculate_log_likelihood( const unsigned long int & position1, const unsigned long int & position2 );
    double calculate_log_likelihood( const vector< double > & );
    double calculate_total_log_likelihood();
    void print_data();
    void get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations); // calculate the number of observations given a set of sufficient statistics

protected:
    //m_data is stored in probability_model (base class)
    bool m_independent;
    double m_markov_alpha;
    unsigned long int m_number_of_states;
};

markov_chain::markov_chain( bool independent, vector< string > filenames, unsigned long int number_of_states, unsigned long int end, unsigned long int diff, double markov_alpha, vector< unsigned long int > separarators):probability_model(end, diff), m_independent(independent), m_markov_alpha(markov_alpha), m_number_of_states(number_of_states) {
    if (m_independent) {
        m_data = new data_input(filenames[0], filenames[1], end, "independent_chain", m_number_of_states, diff);
    } else {
        m_data = new data_input(filenames[0], filenames[1], end, "markov_chain", m_number_of_states, diff, separarators);
    }
}

markov_chain::~markov_chain(){
    m_data->~data_input();
}

double markov_chain::calculate_total_log_likelihood(){
    vector< double > total_cumulative_data;
    m_data->get_cumulative_data_row(m_data->get_data_size() - 1, total_cumulative_data);
    return calculate_log_likelihood(total_cumulative_data);
}

double markov_chain::calculate_log_likelihood(const unsigned long int & position1, const unsigned long int & position2) {
    //return 0;
    vector< double > sufficient_stats1;
    m_data->get_cumulative_data_row(position1, sufficient_stats1);
    vector< double > sufficient_stats2;
    m_data->get_cumulative_data_row(position2, sufficient_stats2);
    for (size_t i = 0; i < sufficient_stats1.size(); i++) {
        sufficient_stats2[i] -= sufficient_stats1[i];
    }
    return calculate_log_likelihood(sufficient_stats2);
}

double markov_chain::calculate_log_likelihood(const vector< double > & sufficient_stats) {
    //return 0;
    // check if there are any observations in the sufficient stats - if so, return 0
    bool any_observations = false;
    unsigned int state = 0;
    while (state < m_number_of_states && !any_observations) {
        any_observations = sufficient_stats[state] > 0.1;
        state++;
    }
    if (!any_observations) {
        return 0;
    }
    double log_like = 0;
    if (m_independent) {
        double n_sum = 0;
        for (unsigned long int state = 0; state < m_number_of_states; state++) {
            log_like += gsl_sf_lngamma(sufficient_stats[state] + m_markov_alpha);
            n_sum += sufficient_stats[state];
        }
        log_like += gsl_sf_lngamma(m_number_of_states * m_markov_alpha) - gsl_sf_lngamma(n_sum + m_number_of_states * m_markov_alpha) - m_number_of_states * gsl_sf_lngamma(m_markov_alpha);
    } else {
        for (unsigned long int state1 = 0; state1 < m_number_of_states; state1++) {
            double n_sum = 0;
            for (unsigned long int state2 = 0; state2 < m_number_of_states; state2++) {
                log_like += gsl_sf_lngamma( sufficient_stats[ m_number_of_states * state1 + state2 ] + m_markov_alpha );
                n_sum += sufficient_stats[ m_number_of_states * state2 + state1 ];
            }
            log_like += gsl_sf_lngamma( m_number_of_states * m_markov_alpha ) - gsl_sf_lngamma( n_sum + m_number_of_states * m_markov_alpha ) - m_number_of_states * gsl_sf_lngamma( m_markov_alpha );
        }
    }
    return log_like;
}

void markov_chain::print_data(){
    m_data->print_data();
}

void markov_chain::get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations) {
    if (sufficient_stats.size() != m_number_of_states) {
        cerr << "sufficient_statistics size doesn't match m_number_of_states" << endl;
    }
    double temp_sum = 0;
    for (unsigned int i = 0; i < m_number_of_states; i++) {
        temp_sum += sufficient_stats[i];
    }
    number_of_observations = temp_sum;
}

#endif
