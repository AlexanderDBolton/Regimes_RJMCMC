#ifndef BINARY_H
#define BINARY_H

class binary {
    
public:
    binary(const double & log_likelihood = 0, const int & left_index = 0);
    double get_log_likelihood() {return m_log_likelihood;}
    void set_log_likelihood(const double & log_likelihood) {m_log_likelihood = log_likelihood;}
    void increase_log_likelihood(const double & increase) {m_log_likelihood += increase;}
    void set_left_index(const int & index) {m_left_index = index;}
    void increase_left_index(const int & change) {m_left_index += change;}
    int get_left_index() {return m_left_index;}
    //void set_right_index(const int & index) {m_right_index = index;}
    //int get_right_index() {return m_right_index;}
    
protected:
    double m_log_likelihood; // the log_likelihood for the interval defined by the binary
    int m_left_index; // the index of the changepoint that bounds the binary to the left
    
};

binary::binary(const double & log_likelihood, const int & left_index):m_log_likelihood(log_likelihood), m_left_index(left_index){
    
}

#endif
