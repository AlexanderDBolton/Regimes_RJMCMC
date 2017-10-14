#ifndef CHANGEPOINT_H
#define CHANGEPOINT_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <map>


class changepoint {
    
    friend bool operator<(const changepoint &, const changepoint &); //compares the positions of changepoints
    friend bool operator>(const changepoint &, const changepoint &);

public:
    changepoint(const unsigned long int & position = 0);
    unsigned long int get_position() const { return m_position; }
    void set_position(const unsigned long int & position) {m_position = position;}
    double get_log_likelihood() const { return m_log_likelihood; }
    void set_log_likelihood( double log_likelihood );
    void set_binary_index(const vector< unsigned int > index) {m_binary_index = index;}
    void set_binary_index_row(const unsigned int process, const unsigned int & index) {m_binary_index[process] = index;}
    vector< unsigned int > get_binary_index() {return m_binary_index;}
    int get_binary_index_row(const size_t & process) {return m_binary_index[process];}
    void increase_binary_index_row(const unsigned int & process, const int & diff) {m_binary_index[process] += diff;}
    void set_full_index(const vector< unsigned int > index) {m_full_index = index;} // needed? Can probably delete
    void set_full_index_row(const unsigned int process, const int & index) {m_full_index[process] = index;}
    vector< unsigned int > get_full_index() const {return m_full_index;}
    int get_full_index_row(const size_t & process) {return m_full_index[process];}
    void increase_full_index_row_for_inserting_unobserved_regime_if_necessary(const unsigned int & process, const unsigned int & unobserved_regime);
    void decrease_full_index_row_for_removing_unobserved_regime_if_necessary(const unsigned int & process, const unsigned int & removed_regime);
    //unsigned int get_regime_index_row(const unsigned int & process) {return m_regime_index[process];}
    void set_full_index_equal_to_binary_index() { m_full_index = m_binary_index; }
    
protected:
    unsigned long int m_position; // the position of the changepoint. If diff > 1 then this will correspond to the position in the contracted time (e.g. diff = 5, actual position = 6 means postion = 1)
    double m_log_likelihood; // the log likelihood in the interval [this changepoint, changepoint to the right)
    vector< unsigned int > m_binary_index; // used for binary marked vectors. The index of the binary that contains this changepoint
    vector< unsigned int > m_full_index; // used for regimes. The index of the regime that contains this changepoint
    //vector< unsigned int > m_regime_index; // stores the regime index for each process for this process
    
};

changepoint::changepoint(const unsigned long int & position):m_position(position) {
    
}

bool operator < ( const changepoint & cp1, const changepoint & cp2 ){
    return ( cp1.get_position() < cp2.get_position() );
}

bool operator > ( const changepoint & cp1, const changepoint & cp2 ){
    return ( cp1.get_position() > cp2.get_position() );
}

void changepoint::set_log_likelihood( double log_likelihood ){
    m_log_likelihood = log_likelihood;
}

void changepoint::increase_full_index_row_for_inserting_unobserved_regime_if_necessary(const unsigned int & process, const unsigned int & unobserved_regime) {
    if (m_full_index[process] >= unobserved_regime) {
        m_full_index[process]++;
    }
}

void changepoint::decrease_full_index_row_for_removing_unobserved_regime_if_necessary(const unsigned int & process, const unsigned int & removed_regime) {
    if (m_full_index[process] > removed_regime) {
        m_full_index[process]--;
    }
}

#endif