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

class mult_process : public probability_model{

  public:

    mult_process(const string, const string, const unsigned long int &, const unsigned long int &, double &);
    ~mult_process();
    double cost(const unsigned long int & position_0, const unsigned long int & position_1);
    void print_data();

  private:

    std::size_t m_number_of_processes; //number of processes
    vector< probability_model * > m_processes;
};

mult_process::mult_process(const string file_of_filenames, const string penalty, const unsigned long int & end_time, const unsigned long int & diff, double & beta):probability_model(end_time, diff) {
    string model_type = "gaussian"; //default model type
    string line, filename;
    ifstream DataFilenames(file_of_filenames, ios::in);
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
            if (model_type == "gaussian") {
                m_processes.push_back(new gaussian_model(filenames, m_end, m_diff));
                if (penalty == "BIC") {
                    beta += 2.0 * log(static_cast< double >(end_time + 1));
                }
                else if (penalty == "AIC") {
                    beta += 2.0 * 2.0;
                }
            }
            if (model_type == "poisson") {
                m_processes.push_back(new poisson_model(filenames, m_end, m_diff));
                if (penalty == "BIC") {
                    beta += 1.0 * log(static_cast< double >(end_time + 1));
                }
                else if (penalty == "AIC") {
                    beta += 1.0 * 2.0;
                }
            }
            if (model_type == "markov_chain") {
                m_processes.push_back(new markov_chain(false, filenames, number_of_states, m_end, m_diff));
                if (penalty == "BIC") {
                    beta += static_cast< double >(number_of_states * (number_of_states - 1)) * log(static_cast< double >(end_time));
                }
                else if (penalty == "AIC") {
                    beta += static_cast< double >(number_of_states * (number_of_states - 1)) * 2.0;
                }
            }
            if (model_type == "independent_chain") {
                m_processes.push_back(new markov_chain(true, filenames, number_of_states, m_end, m_diff));
                if (penalty == "BIC") {
                    beta += static_cast< double >(number_of_states - 1) * log(static_cast< double >(end_time + 1));
                }
                else if (penalty == "AIC") {
                    beta += static_cast< double >(number_of_states - 1) * 2.0;
                }
                
            }
        }
    }
    m_number_of_processes = m_processes.size();
}

mult_process::~mult_process() {
    for (unsigned int i = 0; i < m_processes.size(); i++) {
        if(m_processes[i]!=nullptr)
            delete m_processes[i];
    }
}

double mult_process::cost(const unsigned long int & position_0, const unsigned long int & position_1) {
    double cost = 0;
    for (unsigned int process = 0; process < m_number_of_processes; process++) {
        cost += m_processes[process]->cost(position_0, position_1);
    }
    return cost;
}

void mult_process::print_data() {
    
}

#endif
