#ifndef DATA_INPUT_H
#define DATA_INPUT_H

#include <fstream>
#include <string>
#include <vector>
using namespace std;

class data_input {

public:
    data_input(string data_filename, string data_times_filename = "uniform_arrivals", unsigned long int end_time = 1, string model_type = "gaussian", unsigned long int number_of_states = 3, unsigned long int diff = 1, const vector< unsigned long int > & separators = vector< unsigned long int >(0));
    ~data_input();
    void print_data();
    size_t get_data_size() { return m_cumulative_data.size(); }
    void get_cumulative_data_row(const unsigned long int & index, vector< double > & cumulative_data);

protected:
    //vector< vector< double > > m_cumulative_observations; //for Gaussian data, first column gives cumulative observations, second column gives cumulative number of observations
    vector< vector< double > > m_cumulative_data; //for Markov chains, independent chains and Poisson data. See manual for how the data is fed in.
    unsigned long int m_diff = 1; //the data is compacted so that each record in m_cumulative_observations or m_cumulative_data represents m_diff observations. Makes it easier to perform inference, e.g. only look for changepoints every 10 observations.
    unsigned long int m_number_of_states;
    vector< unsigned long int > m_separators;
};

data_input::data_input(string data_filename, string data_times_filename, unsigned long int end_time, string model_type, unsigned long int number_of_states, unsigned long int diff, const vector< unsigned long int > & separators):m_diff(diff), m_number_of_states(number_of_states), m_separators(separators) {
    if(model_type == "gaussian") {
        /*for (unsigned long int index = 0; index < (end_time + 1) / m_diff; index++) {
            m_cumulative_data.push_back(vector<double>(m_number_of_states, 0.0));
        }*/
        m_cumulative_data = vector< vector< double > >((end_time + 1) / m_diff, vector< double >(m_number_of_states, 0.0));
        if ((end_time + 1) % m_diff != 0) { //if m_diff doesn't evenly divide end_time + 1 then we will need an extra row for the overspill
            m_cumulative_data.push_back(vector< double >(m_number_of_states, 0.0));
        }
        ifstream dataFile(data_filename.c_str(), ifstream::in);
        unsigned long int actual_obs = 0, record_index = 0; //actual_obs keeps track of the actual observation that we are feeding in, obs_index keeps track of where this observation is stored in m_cumulative_data
        double cumulative_sum = 0, cumulative_sum_of_squares = 0, number_of_observations = 0;
        double data_value;
        if (data_times_filename == "uniform_arrivals") {
            while (dataFile >> data_value) {
                cumulative_sum += data_value;
                cumulative_sum_of_squares += data_value * data_value;
                number_of_observations++;
                if (actual_obs % m_diff == m_diff - 1) {
                    vector< double > temp_vec(3, 0.0);
                    temp_vec[0] = cumulative_sum;
                    temp_vec[1] = cumulative_sum_of_squares;
                    temp_vec[2] = number_of_observations;
                    m_cumulative_data[record_index] = temp_vec;
                    cumulative_sum = 0;
                    cumulative_sum_of_squares = 0;
                    number_of_observations = 0;
                    record_index++;
                }
                actual_obs++;
            }
            if (number_of_observations > 0.5) {//if there are some values that haven't been put in M-cumulative_observations because M-diff doesn't evenly divide end_time + 1.
                vector< double > temp_vec(3, 0.0);
                temp_vec[0] = cumulative_sum;
                temp_vec[1] = cumulative_sum_of_squares;
                temp_vec[2] = number_of_observations;
                m_cumulative_data[(end_time + 1) / m_diff] = temp_vec;
            }
        }
        else { //non-uniform arrivals
            ifstream timesFile(data_times_filename.c_str(), ifstream::in);
            unsigned long int time_value, max_time = m_diff;
            while (dataFile >> data_value && timesFile >> time_value) {
                if (time_value >= max_time) {
                    vector< double > temp_vec(3, 0.0);
                    temp_vec[0] = cumulative_sum;
                    temp_vec[1] = cumulative_sum_of_squares;
                    temp_vec[2] = number_of_observations;
                    m_cumulative_data[record_index] = temp_vec;
                    cumulative_sum = 0;
                    cumulative_sum_of_squares = 0;
                    number_of_observations = 0;
                    while (time_value >= max_time) {
                        record_index++;
                        max_time += m_diff;
                    }
                }
                cumulative_sum += data_value;
                cumulative_sum_of_squares += data_value * data_value;
                number_of_observations++;
            }
            if (number_of_observations > 0.5) {//if there are some values that haven't been put in M-cumulative_observations because M-diff doesn't evenly divide end_time + 1.
                vector< double > temp_vec(3, 0.0);
                temp_vec[0] = cumulative_sum;
                temp_vec[1] = cumulative_sum_of_squares;
                temp_vec[2] = number_of_observations;
                m_cumulative_data[record_index] = temp_vec;
            }
        }
    }
    else if (model_type == "poisson") {
        /*for (unsigned long int index = 0; index < (end_time+1)/m_diff; index++) {
            m_cumulative_data.push_back(vector< double >(m_number_of_states, 0.0));
        }*/
        m_cumulative_data = vector< vector< double > >((end_time + 1) / m_diff, vector< double >(m_number_of_states, 0.0));
        if ((end_time + 1) % m_diff != 0) { //if m_diff doesn't evenly divide end_time + 1 then we will need an extra row for the overspill
            m_cumulative_data.push_back(vector< double >(m_number_of_states, 0.0));
        }
        ifstream dataFile(data_filename.c_str(), ifstream::in);
        unsigned long int actual_obs = 0, record_index = 0; //actual_obs keeps track of the actual observation that we are feeding in, obs_index keeps track of where this observation is stored in m_cumulative_observations
        double cumulative_sum = 0, number_of_observations = 0, data_value;
        if (data_times_filename == "uniform_arrivals") {
            while (dataFile >> data_value) {
                cumulative_sum += data_value;
                number_of_observations++;
                if (actual_obs % m_diff == m_diff - 1) {
                    vector<double> temp_vec(2, 0);
                    temp_vec[0] = cumulative_sum;
                    temp_vec[1] = number_of_observations;
                    m_cumulative_data[record_index] = temp_vec;
                    cumulative_sum = 0;
                    number_of_observations = 0;
                    record_index++;
                }
                actual_obs++;
            }
            if (number_of_observations > 0.5) {//if there are some values that haven't been put in M-cumulative_observations because M-diff doesn't evenly divide end_time + 1.
                vector< double > temp_vec(2, 0);
                temp_vec[0] = cumulative_sum;
                temp_vec[1] = number_of_observations;
                m_cumulative_data[(end_time + 1) / m_diff] = temp_vec;
            }
        }
        else { //non-uniform arrivals
            ifstream timesFile(data_times_filename.c_str(), ifstream::in);
            unsigned long int time_value, max_time = m_diff;
            while (dataFile >> data_value && timesFile >> time_value) {
                if (time_value >= max_time) {
                    vector<double> temp_vec(2, 0);
                    temp_vec[0] = cumulative_sum;
                    temp_vec[1] = number_of_observations;
                    m_cumulative_data[record_index] = temp_vec;
                    cumulative_sum = 0;
                    number_of_observations = 0;
                    while (time_value >= max_time) {
                        record_index++;
                        max_time += m_diff;
                    }
                }
                cumulative_sum += data_value;
                number_of_observations++;
            }
            if (number_of_observations > 0.5) {//if there are some values that haven't been assigned
                vector<double> temp_vec(2, 0);
                temp_vec[0] = cumulative_sum;
                temp_vec[1] = number_of_observations;
                m_cumulative_data[record_index] = temp_vec;
            }
        }
    }
    else if (model_type == "markov_chain") {
        /*for (unsigned long int index = 0; index < (end_time + 1) / m_diff; index++) {
            m_cumulative_data.push_back(vector< double >(m_number_of_states * m_number_of_states, 0.0));
        }*/
        m_cumulative_data = vector< vector< double > >((end_time + 1) / m_diff, vector< double >(m_number_of_states, 0.0));
        if ((end_time + 1) % m_diff != 0) { //if m_diff doesn't evenly divide end_time + 1 then we will need an extra row for the overspill
            m_cumulative_data.push_back(vector< double >(m_number_of_states * m_number_of_states, 0.0));
        }
        ifstream dataFile(data_filename.c_str(), ifstream::in);
        unsigned long int actual_obs = 0, record_index = 0; //actual_obs keeps track of the actual observation that we are feeding in, obs_index keeps track of where this observation is stored in m_cumulative_observations
        unsigned long int previous_data_value, data_value;
        unsigned int trace_index = 0;
        vector< double > transitions_vector(m_number_of_states * m_number_of_states, 0);
        if (data_times_filename == "uniform_arrivals") {
            dataFile >> previous_data_value;
            actual_obs++;
            if (m_diff == 1) {
                record_index++;
            }
            while (dataFile >> data_value) {
                if (trace_index < m_separators.size() && actual_obs == m_separators[trace_index]) {
                    trace_index++;
                }
                else {
                    transitions_vector[previous_data_value * m_number_of_states + data_value]++;
                }
                if (actual_obs % m_diff == m_diff - 1) {
                    m_cumulative_data[record_index] = transitions_vector;
                    fill(transitions_vector.begin(), transitions_vector.end(), 0);
                    record_index++;
                }
                actual_obs++;
                previous_data_value = data_value;
            }
            bool any_non_zero_elements_in_transitions_vector = false;
            unsigned long int i = 0;
            while (!any_non_zero_elements_in_transitions_vector & (i < transitions_vector.size())) {
                any_non_zero_elements_in_transitions_vector = (transitions_vector[i] > 0.5);
                i++;
            }
            if (any_non_zero_elements_in_transitions_vector) {//if there are some values that haven't been put in m_cumulative_observations because m_diff doesn't evenly divide end_time + 1.
                m_cumulative_data[record_index] = transitions_vector;
            }
        }
        else { //non-uniform arrivals
            ifstream timesFile(data_times_filename.c_str(), ifstream::in);
            unsigned long int time_value, max_time = m_diff;
            dataFile >> previous_data_value;
            timesFile >> time_value; //intentional repetition since we are not interested in the first time value
            while (dataFile >> data_value && timesFile >> time_value) {
                if (time_value >= max_time) {
                    m_cumulative_data[record_index] = transitions_vector;
                    fill(transitions_vector.begin(), transitions_vector.end(), 0);
                    while (time_value >= max_time){
                        record_index++;
                        max_time += m_diff;
                    }
                }
                transitions_vector[previous_data_value * m_number_of_states + data_value]++;
                previous_data_value = data_value;
            }
            bool any_non_zero_elements_in_transitions_vector = false;
            unsigned long int i = 0;
            while (!any_non_zero_elements_in_transitions_vector & (i < transitions_vector.size())) {
                any_non_zero_elements_in_transitions_vector = (transitions_vector[i] > 0.5);
                i++;
            }
            if (any_non_zero_elements_in_transitions_vector) {//if there are some values that haven't been put in m_cumulative_observations because m_diff doesn't evenly divide end_time + 1.
                m_cumulative_data[record_index] = transitions_vector;
            }
        }
    }
    else if (model_type == "independent_chain") {
        cout << "reading in data" << endl;
        /*for (unsigned long int index = 0; index < (end_time + 1) / m_diff; index++) {
            m_cumulative_data.push_back(vector< double >(m_number_of_states, 0.0));
        }*/
        m_cumulative_data = vector< vector< double > >((end_time + 1) / m_diff, vector< double >(m_number_of_states, 0.0));
        if ((end_time + 1) % m_diff != 0) { //if m_diff doesn't evenly divide end_time + 1 then we will need an extra row for the overspill
            m_cumulative_data.push_back(vector< double >(m_number_of_states, 0.0));
        }
        cout << "finished creating m_cumulative_data" << endl;
        ifstream dataFile(data_filename.c_str(), ifstream::in);
        unsigned long int actual_obs = 0, record_index = 0; //actual_obs keeps track of the actual observation that we are feeding in, obs_index keeps track of where this observation is stored in m_cumulative_observations
        unsigned long int data_value;
        vector< double > observations_vector(m_number_of_states, 0);
        if (data_times_filename == "uniform_arrivals") {
            while (dataFile >> data_value){
                observations_vector[data_value]++;
                if (actual_obs % m_diff == m_diff - 1) {
                    m_cumulative_data[record_index] = observations_vector;
                    fill(observations_vector.begin(), observations_vector.end(), 0);
                    record_index++;
                }
                actual_obs++;
            }
            bool any_non_zero_elements_in_observations_vector = false;
            unsigned long int i = 0;
            while (!any_non_zero_elements_in_observations_vector & (i < observations_vector.size())) {
                any_non_zero_elements_in_observations_vector = (observations_vector[i] > 0.5);
                i++;
            }
            if (any_non_zero_elements_in_observations_vector) {//if there are some values that haven't been put in m_cumulative_observations because m_diff doesn't evenly divide end_time + 1.
                m_cumulative_data[record_index] = observations_vector;
            }
        }
        else { //non-uniform arrivals
            ifstream timesFile(data_times_filename.c_str(), ifstream::in);
            unsigned long int time_value, max_time = m_diff;
            while (dataFile >> data_value && timesFile >> time_value) {
                if (time_value >= max_time) {
                    m_cumulative_data[record_index] = observations_vector;
                    fill(observations_vector.begin(), observations_vector.end(), 0);
                    while (time_value >= max_time) {
                        record_index++;
                        max_time += m_diff;
                    }
                }
                observations_vector[data_value]++;
            }
            bool any_non_zero_elements_in_observations_vector = false;
            unsigned long int i = 0;
            while (!any_non_zero_elements_in_observations_vector & (i < observations_vector.size())) {
                any_non_zero_elements_in_observations_vector = (observations_vector[i] > 0.5);
                i++;
            }
            if (any_non_zero_elements_in_observations_vector) {//if there are some values that haven't been put in m_cumulative_observations because m_diff doesn't evenly divide end_time + 1.
                m_cumulative_data[record_index] = observations_vector;
            }
        }
    }
    else {
        cout << "no model type was provided" << '\n';
    }
    
    //insert vector of 0s at the beginning
    m_cumulative_data.insert(m_cumulative_data.begin(), vector< double >(m_number_of_states, 0.0));

    //convert data into cumulative data
    size_t cumulative_data_size = m_cumulative_data.size();
    for (size_t i = 1; i < cumulative_data_size; i++) {
        for (size_t j = 0; j < m_number_of_states; j++) {
            m_cumulative_data[i][j] += m_cumulative_data[i - 1][j];
        }
    }
}

data_input::~data_input(){
    
}

void data_input::print_data(){
  for (unsigned long int i = 0; i < m_cumulative_data.size(); i++) {
    for (int j = 0; j < m_cumulative_data[i].size(); j++) {
      cout << m_cumulative_data[i][j] << " ";
    }
    cout << '\n';
  }
}

void data_input::get_cumulative_data_row(const unsigned long int & index, vector< double > & cumulative_data) {
    if (index < m_cumulative_data.size()) {
        cumulative_data = m_cumulative_data[index];
    } else {
        cumulative_data = m_cumulative_data[m_cumulative_data.size() - 1];
    }
}


#endif
