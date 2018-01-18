#ifndef MARKOV_CHAIN_MODEL_H
#define MARKOV_CHAIN_MODEL_H

#include <math.h>
#include <stdlib.h>
using namespace std;

class markov_chain : public probability_model{

public:
    markov_chain(const bool &, const vector< string > &, const unsigned int &, const unsigned long int &, const unsigned long int &);
    ~markov_chain() {}
    void get_cumulative_sufficient_data(const unsigned long int & time, vector< double > & sufficient_data);
    double cost(const unsigned long int &, const unsigned long int &); // calculate minus the log likelihood for observations betwixt position_0 and position_1 (inclusive).
    vector< double > calculate_sufficient_statistics(const unsigned long int &, const unsigned long int &); // calculate the sufficient statistics for an interval, including both end points
    vector< double > calculate_maximum_likelihood_probs(const vector< double > &, const unsigned int & = 0); // the offset parameter is used when the process is not independent. As the data is stored in the format (n_00, n_01, ..., n_0c, n_10, n_11, ..., n_1c, ...) we need to know which state we are on so that we can start going through the vector at element offset * number_of_states.
    void print_data();

protected:
    bool m_independent;
    unsigned int m_number_of_states;
    data_input m_data;
};

markov_chain::markov_chain(const bool & independent, const vector< string > & filenames, const unsigned int & number_of_states, const unsigned long int & end, const unsigned long int & diff):probability_model(end, diff), m_independent(independent), m_number_of_states(number_of_states), m_data(filenames[0], filenames[1], end, independent ? "independent_chain" : "markov_chain", m_number_of_states, diff) {
}

void markov_chain::get_cumulative_sufficient_data(const unsigned long int & time, vector< double > & sufficient_data){
    m_data.get_cumulative_data_row(time, sufficient_data);
}

double markov_chain::cost(const unsigned long int & position_0, const unsigned long int & position_1) {
    vector< double > interval_sufficient_statistics = calculate_sufficient_statistics(position_0, position_1);
    
    // calculate the log likelihood of the observations
    double maximum_log_likelihood = 0;
    if (m_independent) {
        vector< double > maximum_likelihood_probs = calculate_maximum_likelihood_probs(interval_sufficient_statistics);
        for (unsigned int i = 0; i < m_number_of_states; i++) {
            if (0.5 < interval_sufficient_statistics[i]) {
                maximum_log_likelihood += interval_sufficient_statistics[i] * log(maximum_likelihood_probs[i]);
            }
        }
    }
    else {
        for (unsigned int i = 0; i < m_number_of_states; i++) {
            vector< double > maximum_likelihood_probs = calculate_maximum_likelihood_probs(interval_sufficient_statistics, i);
            for (unsigned int j = 0; j < m_number_of_states; j++) {
                if (0.5 < interval_sufficient_statistics[m_number_of_states * i + j]) {
                    maximum_log_likelihood += interval_sufficient_statistics[m_number_of_states * i + j] * log(maximum_likelihood_probs[j]);
                }
            }
        }
    }
    return -maximum_log_likelihood * 2.0;
}

vector< double > markov_chain::calculate_sufficient_statistics(const unsigned long int & position_0, const unsigned long int & position_1) {
    // calculate sufficient statistics for this interval
    vector< double > sufficient_stats_0;
    m_data.get_cumulative_data_row(position_0, sufficient_stats_0);
    vector< double > sufficient_stats_1;
    m_data.get_cumulative_data_row(position_1, sufficient_stats_1);
    for (unsigned int i = 0; i < sufficient_stats_1.size(); i++) {
        sufficient_stats_1[i] -= sufficient_stats_0[i];
    }
    return sufficient_stats_1;
}

vector< double > markov_chain::calculate_maximum_likelihood_probs(const vector< double > & sufficient_statistics, const unsigned int & offset) {
    // calculate the sum of the sufficient statistics
    double sum = 0;
    for (unsigned int index = offset * m_number_of_states; index < offset * m_number_of_states + m_number_of_states; index++) {
        sum += sufficient_statistics[index];
    }
    // now calculate the vector
    vector< double > ml_probabilities = vector< double >(m_number_of_states);
    for (unsigned int index = offset * m_number_of_states; index < offset * m_number_of_states + m_number_of_states; index++) {
        ml_probabilities[index - offset * m_number_of_states] = sufficient_statistics[index] / sum;
    }
    return ml_probabilities;
}

void markov_chain::print_data() {
    m_data.print_data();
}

#endif
