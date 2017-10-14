#ifndef PROBABILITY_MODEL_H
#define PROBABILITY_MODEL_H

#include <vector>
#include <gsl/gsl_sf_gamma.h>
#include "data_input.h"
using namespace std;

class probability_model {

public:
	probability_model( unsigned long int end, unsigned long int diff);
    virtual ~probability_model();
    virtual double calculate_log_likelihood( const unsigned long int & position1, const unsigned long int & position2 ) = 0;
    virtual double calculate_log_likelihood( const vector< double > & ) = 0;
    unsigned long int get_diff() {return m_diff;}
    virtual void print_data() = 0;
    //virtual double calculate_total_log_likelihood() = 0;
    void get_cumulative_sufficient_data(const unsigned long int & time, vector< double > & sufficient_data);
    virtual void get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations) = 0; // calculate the number of observations given a set of sufficient statistics
    unsigned long int m_end;
    unsigned long int m_diff;
    data_input * m_data;
};

probability_model::probability_model( unsigned long int end, unsigned long int diff ):m_end( end ), m_diff( diff){

}

probability_model::~probability_model(){
    
}

void probability_model::get_cumulative_sufficient_data(const unsigned long int & time, vector< double > & sufficient_data){
    m_data->get_cumulative_data_row(time, sufficient_data);
}

#endif