#ifndef GAUSSIAN_MODEL_REGIME_h
#define GAUSSIAN_MODEL_REGIME_h

class gaussian_model: public probability_model{
    
public:
    gaussian_model(vector< string > filenames, unsigned long int end, unsigned long int diff);
    ~gaussian_model() {}
    void get_cumulative_sufficient_data(const unsigned long int & time, vector< double > & sufficient_data);
    double cost(const unsigned long int &, const unsigned long int &);
    vector< double > calculate_sufficient_statistics(const unsigned long int &, const unsigned long int &);
    void print_data();
    
protected:
    double m_mu;
    double m_tau;
    data_input m_data;
};

gaussian_model::gaussian_model(vector< string > filenames, unsigned long int end, unsigned long int diff):probability_model(end, diff), m_data(filenames[0], filenames[1], end, "gaussian", 3, diff) {
    //m_data = new data_input(filenames[0], filenames[1], end, "gaussian", 3, diff);
}

void gaussian_model::get_cumulative_sufficient_data(const unsigned long int & time, vector< double > & sufficient_data) {
    m_data.get_cumulative_data_row(time, sufficient_data);
}

double gaussian_model::cost(const unsigned long int & position_0, const unsigned long int & position_1) {
    vector< double > interval_sufficient_statistics = calculate_sufficient_statistics(position_0, position_1);
    double sum_x = interval_sufficient_statistics[0];
    double sum_x_squared = interval_sufficient_statistics[1];
    double n = interval_sufficient_statistics[2];
    
    // check if there is only one observation in the interval
    if (n < 1.5) {
        // the likelihood for the observation will be 1, so log(likelihood) = 0
        return 0;
    }
    
    double mu_hat = sum_x / n;
    double sigma_squared_hat = (sum_x_squared / n) - (mu_hat * mu_hat);
    
    if (sigma_squared_hat < 0) {
        cerr << "Gaussian variance estimated to be below 0";
    }
    
    // check if the variance is 0
    if (sigma_squared_hat < 0.0000001) {
        // the likelihood for the observation will be 1, so log(likelihood) = 0
        return 0;
    }
    
    double maximum_log_likelihood = -(log(2*M_PI) + log(sigma_squared_hat) + 1) * n / 2;
    return -maximum_log_likelihood * 2.0;
}

vector< double > gaussian_model::calculate_sufficient_statistics(const unsigned long int & position_0, const unsigned long int & position_1) {
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



void gaussian_model::print_data() {
    m_data.print_data();
}

#endif
