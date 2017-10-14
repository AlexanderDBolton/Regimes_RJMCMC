#ifndef MULTIPLE_PROCESSES_REGIME_H
#define MULTIPLE_PROCESSES_REGIME_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "basic_probability_model.h"
//Insert any probability models to be used here:
#include "poisson_model_regime.h"
#include "gaussian_model_regime.h"
#include "markov_chain_model_regime.h"
//#include "changepoint.h"

class mult_process : public probability_model{

  public:

    mult_process(const unsigned long int &, const unsigned long int &, const string, const vector< double > &, const vector< unsigned long int > &);
    ~mult_process();
    double calculate_total_log_likelihood();
    double calculate_log_likelihood( const unsigned int & process, const vector<double> & sufficient_statistics );
	double calculate_log_likelihood(const vector< vector< double > > & sufficient_statistics);
    double calculate_log_likelihood( const unsigned long int & position1, const unsigned long int & position2 );
    double calculate_log_likelihood(const unsigned int & process, const unsigned long int & position1, const unsigned long int & position2);
    void get_cumulative_sufficient_data(const unsigned long int & index, vector< vector< double > > & sufficient_stats); // sufficient_stats will be assigned to be the vector of vectors of sufficient statistics as time
    void get_number_of_observations(const vector< vector< double > > & sufficient_stats, vector< double > & number_of_observations);
    virtual void get_number_of_observations(const vector< double > & sufficient_stats, double & number_of_observations){}
    void get_number_of_observations(const unsigned int process, const vector< double > & sufficient_statistics, double & number_of_observations);
    std::size_t get_number_of_processes() { return m_number_of_processes; }
  //virtual double log_likelihood_changepoints( unsigned int, vector<unsigned long int>&, vector<double>& );
    virtual double calculate_log_likelihood(const vector< double > & sufficient_statistics) { return 0; }
    void print_data();

  private:

    std::size_t m_number_of_processes; //number of processes
    vector< probability_model* > m_processes;
    //double m_likelihood_term;
    //string m_model_type;
};

mult_process::mult_process(const unsigned long int & end_time, const unsigned long int & diff, const string file_of_filenames, const vector< double > & probability_model_parameters, const vector< unsigned long int > & separators):probability_model(end_time, diff) {
    string model_type = "gaussian"; //default model type
    string line, filename;
    ifstream DataFilenames( file_of_filenames, ios::in );
    while (DataFilenames.good()) {
        getline(DataFilenames, line);
        vector< string > filenames;
        if (line.size() > 0) {
            istringstream iss(line);
            unsigned int number_of_states;
            iss >> number_of_states;
            string model_type;
            iss >> model_type;
            while (iss >> filename)
                filenames.push_back(filename);
            if (model_type == "gaussian")
                m_processes.push_back(new gaussian_model(filenames, m_end, m_diff, probability_model_parameters[0], probability_model_parameters[1], probability_model_parameters[2], probability_model_parameters[3]));
            if (model_type == "poisson")
                m_processes.push_back(new poisson_model(filenames, m_end, m_diff, probability_model_parameters[4], probability_model_parameters[5]));
            if (model_type == "markov_chain")
                m_processes.push_back(new markov_chain(false, filenames, number_of_states, m_end, m_diff, probability_model_parameters[6], separators));
            if (model_type == "independent_chain")
                m_processes.push_back(new markov_chain(true, filenames, number_of_states, m_end, m_diff, probability_model_parameters[6]));
        }
    }
    m_number_of_processes = m_processes.size();
}

mult_process::~mult_process() {
    for (std::size_t i = 0; i < m_number_of_processes; i++) {
        delete m_processes[i];
    }
}

/*double mult_process::calculate_total_log_likelihood(){
    double total_log_likelihood = 0;
    for( size_t process = 0; process < m_processes.size(); process++ ){
        total_log_likelihood += m_processes[ process ]->calculate_total_log_likelihood();
    }
    return total_log_likelihood;
}*/
                                      
double mult_process::calculate_log_likelihood(const unsigned int & process, const vector<double> & sufficient_statistics) {
    return m_processes[process]->calculate_log_likelihood(sufficient_statistics);
}

// calculate the log likelihood over all processes given sufficient statistics
double mult_process::calculate_log_likelihood(const vector< vector< double > > & sufficient_statistics) {
	double log_likelihood = 0;
	for (unsigned int process = 0; process < m_number_of_processes; process++) {
		log_likelihood += m_processes[process]->calculate_log_likelihood(sufficient_statistics[process]);
	}
	return log_likelihood;
}

double mult_process::calculate_log_likelihood(const unsigned long int & position1, const unsigned long int & position2) {
    double log_likelihood = 0;
    for (size_t process = 0; process < m_processes.size(); process++) {
        log_likelihood += m_processes[process]->calculate_log_likelihood(position1, position2);
    }
    return log_likelihood;
}

double mult_process::calculate_log_likelihood(const unsigned int & process, const unsigned long int & position1, const unsigned long int & position2) {
    return m_processes[process]->calculate_log_likelihood(position1, position2);
}
                                      
void mult_process::print_data(){
    for (size_t process = 0; process < m_number_of_processes; process++) {
        m_processes[process]->print_data();
    }
}

// the sufficient_stats vector must be initialised as a vector< vector< double > >(m_number_of_processes) before being fed in. The function assigns it to be the vector of vectors of sufficient statistics
void mult_process::get_cumulative_sufficient_data(const unsigned long int & index, vector< vector< double > > & sufficient_stats) {
    for (unsigned int process = 0; process < m_number_of_processes; process++){
        m_processes[process]->get_cumulative_sufficient_data(index, sufficient_stats[process]);
    }
}

void mult_process::get_number_of_observations(const vector< vector< double > > & sufficient_stats, vector< double > & number_of_observations) {
    for (unsigned int process = 0; process < m_number_of_processes; process++) {
        m_processes[process]->get_number_of_observations(sufficient_stats[process], number_of_observations[process]);
    }
}

void mult_process::get_number_of_observations(const unsigned int process, const vector< double > & sufficient_statistics, double & number_of_observations) {
    m_processes[process]->get_number_of_observations(sufficient_statistics, number_of_observations);
}

/*double mult_process::log_likelihood_interval(changepoint *obj1, changepoint *obj2){

    double like=0;
  

    for(unsigned int i=0; i<m_k; i++){
      obj1->set_index_from_int(i);
      obj2->set_index_from_int(i);
      like+=m_k_processes[i]->log_likelihood_interval(obj1,obj2);
    }
    return(like);
}


void mult_process::set_data_index(changepoint * cpobj, unsigned int row, changepoint * cpobj_left, changepoint * cpobj_right)
{
  int * vec = NULL;
  if(m_k)
    vec = new int[m_k];
  for(unsigned int i=0; i<m_k; i++){
    if(cpobj_left)
      cpobj_left->set_index_from_int(i);
    if(cpobj_right)
      cpobj_right->set_index_from_int(i);
    m_k_processes[i]->set_data_index(cpobj,row,cpobj_left,cpobj_right);
    vec[i]=cpobj->getdataindex();
  }
  cpobj->set_int_vector(vec,m_k);
}

double mult_process::log_likelihood_changepoints( unsigned int pro, vector<unsigned long int>& regime_changepoints_data_indices, vector<double>& regime_changepoints_changepoint_positions ){
  //set the data indices for the chosen CPs.
  return m_k_processes[ pro ]->log_likelihood_changepoints( regime_changepoints_data_indices, regime_changepoints_changepoint_positions );
}*/

#endif
