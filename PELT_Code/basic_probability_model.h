#ifndef PROBABILITY_MODEL_H
#define PROBABILITY_MODEL_H

#include <vector>
#include "data_input.h"
using namespace std;

class probability_model {

public:
    virtual ~probability_model() {};
    virtual void print_data() = 0;
    virtual double cost(const unsigned long int &, const unsigned long int &) = 0;
    unsigned long int get_diff() {return m_diff;}
    
protected:
    probability_model(unsigned long int end, unsigned long int diff);
    unsigned long int m_end;
    unsigned long int m_diff;
};

probability_model::probability_model(unsigned long int end, unsigned long int diff):m_end(end), m_diff(diff) {

}

#endif