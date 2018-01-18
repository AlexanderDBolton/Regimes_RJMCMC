#ifndef POISSON_MODEL_REGIME_h
#define POISSON_MODEL_REGIME_h

class poisson_model: public probability_model {

public:
    poisson_model(vector< string > filenames, unsigned long int end, unsigned long int diff); //construct a Poisson model with shape alpha and rate beta, calls a cumulative_data constructor
    ~poisson_model() {} //destructor for poisson_model
    double cost(const unsigned long int &, const unsigned long int &);
    vector< double > calculate_sufficient_statistics(const unsigned long int &, const unsigned long int &);
    void print_data(); //prints the cumulative data for the Poisson sequence.
    
protected:
    data_input m_data;
};

poisson_model::poisson_model(vector< string > filenames, unsigned long int end, unsigned long int diff):probability_model(end, diff), m_data(data_input(filenames[0], filenames[1], end, "poisson", 2, diff)) {
}

double poisson_model::cost(const unsigned long int & position_0, const unsigned long int & position_1) {
    vector< double > interval_sufficient_statistics = calculate_sufficient_statistics(position_0, position_1);
    double sum_x = interval_sufficient_statistics[0];
    double n = interval_sufficient_statistics[1];
    double sum_of_log_x_factorial = interval_sufficient_statistics[2];
    
    double x_bar = sum_x / n;
    
    double maximum_log_likelihood = (log(x_bar) * sum_x) - (x_bar * n) - sum_of_log_x_factorial;
    return -maximum_log_likelihood * 2.0;
}

vector< double > poisson_model::calculate_sufficient_statistics(const unsigned long int & position_0, const unsigned long int & position_1) {
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

void poisson_model::print_data(){
    m_data.print_data();
}

#endif
