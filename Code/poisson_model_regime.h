#ifndef POISSON_MODEL_REGIME_h
#define POISSON_MODEL_REGIME_h

#include "basic_probability_model.h"
#include "data_input.h"

class poisson_model: public probability_model{

public:
    poisson_model( vector< string > filenames, unsigned long int end, unsigned long int diff, double alpha, double beta ); //construct a Poisson model with shape alpha and rate beta, calls a cumulative_data constructor
    ~poisson_model(); //destructor for poisson_model
    double calculate_log_likelihood( const unsigned long int & position1, const unsigned long int & position2 ); //calculate the marginal log likelihood between position1 and position2. Uses cumulative data[position1] and cumulative_data[position2 - 1] as cumulative_data begins with 0s. Calculates log likelihood by getting sufficient statistics and then passing them into log_likelihood function for sufficient statistics.
    double calculate_log_likelihood( const vector< double > & ); //calculate the likelihood given sufficient statistics.
    void print_data(); //prints the cumulative data for the Poisson sequence.
    void get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations); // calculate the number of observations given a set of sufficient statistics
    
protected:
    double m_alpha; //shape parameter
    double m_beta; //rate parameter
    //m_data is stored in probability_model (base class)
};

poisson_model::poisson_model( vector< string > filenames, unsigned long int end, unsigned long int diff, double alpha, double beta ):probability_model( end, diff ), m_alpha( alpha ), m_beta( beta ){
    m_data = new data_input(filenames[ 0 ], filenames[ 1 ], end, "poisson", 2, diff);
}

poisson_model::~poisson_model(){
    m_data->~data_input();
}

double poisson_model::calculate_log_likelihood( const unsigned long int & position1, const unsigned long int & position2 ){
    //calculate sufficient statistics
    vector< double > sufficient_stats1;
    m_data->get_cumulative_data_row(position1, sufficient_stats1);
    vector< double > sufficient_stats2;
    m_data->get_cumulative_data_row(position2, sufficient_stats2);
    for( size_t i = 0; i < sufficient_stats1.size(); i++ ){
        sufficient_stats2[ i ] -= sufficient_stats1[ i ];
    }
    //feed sufficient statistics into the log likelihood function for sufficient statistics
    return calculate_log_likelihood( sufficient_stats2 );
}

double poisson_model::calculate_log_likelihood( const vector< double > & sufficient_stats ){ //sufficient stats will have length 2, Poisson observations, and number of Poisson RVs observed.
    if (sufficient_stats[1] < 0.1) { // if there were no observations
        return 0;
    }
    return m_alpha * log( m_beta ) + gsl_sf_lngamma( m_alpha + sufficient_stats[ 0 ] ) - gsl_sf_lngamma( m_alpha ) - ( m_alpha + sufficient_stats[ 0 ] ) * log( m_beta + sufficient_stats[ 1 ] );
}

void poisson_model::print_data(){
    m_data->print_data();
}

void poisson_model::get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations) {
    number_of_observations = sufficient_stats[1];
}

#endif
