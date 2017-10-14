#ifndef REGIME_H
#define REGIME_H

#include <limits>

//regime needs to contain the sufficient statistics
class regime {
    
public:
    regime(const vector< int > & right_indices, const vector< unsigned int > & right_transitions, const vector< unsigned int > & right_transitions_histogram, const vector< double > & statistic, const size_t & number_of_traces, const unsigned int & trace_index, const double & number_of_obs);
    regime(){}
    vector< double > get_sufficient_statistics(); // returns the vector of sufficient statistics (possibly of length 1)
    void alter_sufficient_statistic(const vector< double > & change ); // increase the sufficient statistic by change
    size_t get_size_of_sufficient_statistics() { return m_sufficient_statistics.size(); }
    // size_t get_number_of_changepoints(); not needed?
    void add_interval(const int & right_cp_position, vector< bool > & unobserved_regimes, const unsigned int & unobserved_regimes_index, unsigned int & number_of_unobserved_regimes, const unsigned int & transition = -1); // add these changepoint indices to the collection of right indices. Also add this transition to m_right_transitions and m_right_transitions_histogram
    void remove_interval(const int & right_cp_position, vector< bool > & unobserved_regimes, const unsigned int & unobserved_regimes_index, unsigned int & number_of_unobserved_regimes); // remove the right cp_index from m_right_changepoint_indices, remove the necessary transition from m_right_transitions and the histogram
    unsigned long int find_cp_position_in_m_right_changepoint_indices(const int & cp_position);
    double get_log_likelihood() { return m_log_likelihood; }
    void set_log_likelihood(const double & log_likelihood) { m_log_likelihood = log_likelihood; }
    void set_transitions_out(const unsigned int & transitions_out) { m_transitions_out = transitions_out; }
    unsigned int get_transitions_out() { return m_transitions_out; }
    void increase_transitions_out(const unsigned int & transitions_out_change) { m_transitions_out += transitions_out_change; }
    vector< int > get_right_changepoint_indices() { return m_right_changepoint_indices; }
    vector< unsigned int > get_right_transitions() { return m_right_transitions; }
    vector< unsigned int > get_right_transitions_histogram() { return m_right_transitions_histogram; }
    unsigned int get_right_transitions_histogram_element(const unsigned int & index) { return m_right_transitions_histogram[index]; }
    size_t get_number_of_left_indices();
    unsigned int get_number_of_right_transitions();
    double get_number_of_observations(unsigned int trace) { return m_number_of_observations_in_traces[trace]; }
    void add_regime_to_transitions_histogram(const unsigned int & new_regime);
    void add_right_transition(const unsigned int & new_regime, const unsigned int & index, const bool & end_of_transitions = false);
    void remove_right_transition(const unsigned int & index, const bool & end_of_transitions = false);
    void alter_right_transition(const unsigned int & cp_index, const unsigned int & regime);
    void alter_last_transition(const unsigned int & new_regime);
    void append_right_changepoint_indices(const vector< int > & right_changepoint_indices);
    void append_right_transitions(const vector< unsigned int > & right_transitions);
    void add_right_transitions_histogram(const vector< unsigned int > & right_transitions_histogram);
    void insert_new_unobserved_regime(const unsigned int & new_unobserved_regime_index);
    void remove_unobserved_regime(const unsigned int & regime_to_remove);
    void add_sufficient_statistics(const vector< double > & sufficient_statistics);
    void remove_sufficient_statistics(const vector< double > & sufficient_statistics);
    void add_observations(const unsigned int & trace_index, const double & number_of_observations);
    void remove_observations(const unsigned int & trace_index, const double & number_of_observations);
    void increase_right_indices_greater_than_or_equal_to(const unsigned int & index);
    void decrease_right_indices_greater_than(const unsigned int & index);
    void add_to_similarity_matrix(vector< vector< double > > & similarity_matrix);
    void add_to_min_proportion_similarity_matrices(vector< vector< double > > & min_proportion_similarity_matrix_0, vector< vector< double > > & min_proportion_similarity_matrix_1, const vector< unsigned long int > & actual_number_of_observations_in_traces);
    void add_new_trace();
	void permute_transitions_histogram(const vector< unsigned int > & permutation_vector);
	void permute_right_transitions(const vector< unsigned int > & permutation_vector);
    
protected:
    vector< int > m_right_changepoint_indices; // vector (possibly of length 0) of indices which are the right indices for this regime
    vector< unsigned int > m_right_transitions; // vector of the regimes which are transitioned to at the right cp indices
    vector< unsigned int > m_right_transitions_histogram; // histogram of the right transitions from this regime. The same length as the number of regimes for the process that the regime affects.
    vector< double > m_sufficient_statistics; // vector of sufficient statistics for this regime
    vector< double > m_number_of_observations_in_traces; // for each trace, how many observations are made from this regime in this trace?
    double m_log_likelihood;
    unsigned int m_transitions_out;
};

regime::regime(const vector< int > & right_indices, const vector< unsigned int > & right_transitions, const vector< unsigned int > & right_transitions_histogram, const vector< double > & statistic, const size_t & number_of_traces, const unsigned int & trace_index, const double & number_of_obs){ // trace_index gives the trace that contains the regime. A new regime may only belong to a single trace because of the separators at the end of each trace.
    m_right_changepoint_indices = right_indices;
    m_right_transitions = right_transitions;
    m_right_transitions_histogram = right_transitions_histogram;
    m_sufficient_statistics = statistic;
    m_number_of_observations_in_traces = vector< double > (number_of_traces, 0);
    m_number_of_observations_in_traces[trace_index] = number_of_obs;
}

vector< double > regime::get_sufficient_statistics(){
    return m_sufficient_statistics;
}

void regime::alter_sufficient_statistic(const vector< double > & change){
    for (unsigned long int i = 0; i < m_sufficient_statistics.size(); i++ ){
        m_sufficient_statistics[i] += change[i];
    }
}

void regime::add_interval(const int & right_cp_position, vector< bool > & unobserved_regimes, const unsigned int & unobserved_regimes_index, unsigned int & number_of_unobserved_regimes, const unsigned int & transition) {
    // find where the place for the left index among the left indices
    unsigned long int right_index_position;
    if (m_right_changepoint_indices.size() == 0) {
        right_index_position = 0;
    }
    else {
        int min = 0, max = static_cast< int >(m_right_changepoint_indices.size() - 1);
        while (max >= min) {
            int mid = (min + max) / 2;
            if (m_right_changepoint_indices[mid] < right_cp_position) {
                min = mid + 1;
            }
            else if (m_right_changepoint_indices[mid] > right_cp_position) {
                max = mid - 1;
            }
            else if (m_right_changepoint_indices[mid] == right_cp_position) {
                cerr << "changepoint already exists in m_right_changepoint_indices" << endl;
            }
        }
        if (min >= m_right_changepoint_indices.size()) {
            min = static_cast< int >(m_right_changepoint_indices.size() - 1);
        }
        if (m_right_changepoint_indices[min] > right_cp_position) {
            right_index_position = min;
        }
        else {
            right_index_position = min + 1;
        }
    }
    // insert the cp into m_right_changepoint_indices
    m_right_changepoint_indices.insert(m_right_changepoint_indices.begin() + right_index_position, right_cp_position);
    // add the transition to m_right_transitions and to the transitions histogram
    m_right_transitions.insert(m_right_transitions.begin() + right_index_position, transition);
    if (transition != numeric_limits< unsigned int >::max()) {
        m_right_transitions_histogram[transition]++;
    }
    // if this regime was previously unobserved then set it to be observed and decrease the number of unobserved regimes
    if (unobserved_regimes[unobserved_regimes_index]) {
        unobserved_regimes[unobserved_regimes_index] = false;
        number_of_unobserved_regimes--;
    }
}

void regime::remove_interval(const int & right_cp_position, vector< bool > & unobserved_regimes, const unsigned int & unobserved_regimes_index, unsigned int & number_of_unobserved_regimes) {
    unsigned long int right_index_position = find_cp_position_in_m_right_changepoint_indices(right_cp_position);
    if (right_index_position >= m_right_transitions.size()) {
        cerr << "no transition to remove" << endl;
    }
    // remove the transition from m_right_transitions and from the transitions histogram
    if (m_right_transitions[right_index_position] != numeric_limits< unsigned int >::max()) {
        m_right_transitions_histogram[m_right_transitions[right_index_position]]--;
    }
    m_right_transitions.erase(m_right_transitions.begin() + right_index_position);
    // remove the cp from m_right_changepoint_indices
    m_right_changepoint_indices.erase(m_right_changepoint_indices.begin() + right_index_position);

    // if removing this interval makes the regime unobserved, set it to be unobserved and increase the number of unobserved regimes
    if (m_right_changepoint_indices.size() == 0) {
        unobserved_regimes[unobserved_regimes_index] = true;
        number_of_unobserved_regimes++;
    }
}

unsigned long int regime::find_cp_position_in_m_right_changepoint_indices(const int & cp_position) {
    unsigned long int min = 0, max = m_right_changepoint_indices.size() - 1;
    while (max >= min) {
        unsigned long int mid = (min + max) / 2;
        if (m_right_changepoint_indices[mid] == cp_position) {
            return mid;
        }
        else if (m_right_changepoint_indices[mid] < cp_position) {
            min = mid + 1;
        }
        else {
            max = mid - 1;
        }
    }
    cerr << "could not find that cp_position in m_right_changepoint_indices" << endl;
    return 0;
}

size_t regime::get_number_of_left_indices() {
    return m_right_changepoint_indices.size();
}

unsigned int regime::get_number_of_right_transitions() {
    unsigned int number_of_right_transitions = 0;
    for (unsigned int i = 0; i < m_right_transitions_histogram.size(); i++) {
        number_of_right_transitions += m_right_transitions_histogram[i];
    }
    return number_of_right_transitions;
}

void regime::add_regime_to_transitions_histogram(const unsigned int & new_regime) {
    m_right_transitions_histogram.insert(m_right_transitions_histogram.begin() + new_regime, 0);
}

void regime::add_right_transition(const unsigned int & new_regime, const unsigned int & index, const bool & end_of_transitions) {
    if (end_of_transitions) {
        m_right_transitions.push_back(new_regime);
    }
    else {
        m_right_transitions.insert(m_right_transitions.begin() + index, new_regime);
    }
    m_right_transitions_histogram[new_regime]++;
}

void regime::remove_right_transition(const unsigned int & index, const bool & end_of_transitions) {
    if (end_of_transitions) {
        m_right_transitions_histogram[m_right_transitions.back()]--;
        m_right_transitions[m_right_transitions.size() - 1] = numeric_limits< unsigned int >::max();
    }
    else {
        unsigned long int min = 0, max = m_right_changepoint_indices.size() - 1;
        unsigned long int right_cp_index = -1;
        while (max >= min) {
            unsigned long int mid = (min + max) / 2;
            if (m_right_changepoint_indices[mid] == index) {
                right_cp_index = mid;
                min = max + 1;
            }
            else if (m_right_changepoint_indices[mid] < index) {
                min = mid + 1;
            }
            else {
                max = mid - 1;
            }
        }
        
        // was any element of m_right_changepoint_indices equal to index?
        if (right_cp_index == numeric_limits< unsigned int >::max()) {
            cerr << "no right changepoints found to match index" << endl;
        }
        else {
            // now alter the corresponding element of m_right_transitions
            if (right_cp_index >= m_right_transitions.size()) {
                cerr << "no right transition to remove" << endl;
            }
            m_right_transitions_histogram[m_right_transitions[right_cp_index]]--;
            m_right_transitions[right_cp_index] = numeric_limits< unsigned int >::max();
        }
    }
}

void regime::alter_right_transition(const unsigned int & cp_index, const unsigned int & regime) {
    // find the element of the right changepoint indices that is equal to index
    unsigned long int min = 0, max = m_right_changepoint_indices.size() - 1;
    unsigned long int right_cp_index = -1;
    bool found = false;
    while (max >= min && !found) {
        unsigned long int mid = (min + max) / 2;
        if (m_right_changepoint_indices[mid] == cp_index) {
            right_cp_index = mid;
            found = true;
        }
        else if (m_right_changepoint_indices[mid] < cp_index) {
            min = mid + 1;
        }
        else {
            max = mid - 1;
        }
    }
    
    // was any element of m_right_changepoint_indices equal to index?
    if (right_cp_index == numeric_limits< unsigned long int >::max()) {
        cerr << "no right changepoints found to match index" << endl;
    }
    else {
        // now alter the corresponding element of m_right_transitions
        if (right_cp_index >= m_right_transitions.size()) {
            cerr << "no right transition to replace" << endl;
        }
        if (m_right_transitions[right_cp_index] != numeric_limits< unsigned int >::max()) {
            m_right_transitions_histogram[m_right_transitions[right_cp_index]]--;
        }
        else {
            m_transitions_out++;
        }
        m_right_transitions[right_cp_index] = regime;
        m_right_transitions_histogram[regime]++;
    }
}

void regime::alter_last_transition(const unsigned int & new_regime) {
    unsigned int last_transition = m_right_transitions.back();
    if (last_transition != numeric_limits< unsigned int >::max()) {
        m_right_transitions_histogram[last_transition]--;
    }
    else {
        m_transitions_out++;
    }
    m_right_transitions.back() = new_regime;
    m_right_transitions_histogram[new_regime]++;
}

void regime::append_right_changepoint_indices(const vector< int > & right_changepoint_indices) {
    for (unsigned int i = 0; i < right_changepoint_indices.size(); i++) {
        m_right_changepoint_indices.push_back(right_changepoint_indices[i]);
    }
}

void regime::append_right_transitions(const vector< unsigned int > & right_transitions) {
    for (unsigned int i = 0; i < right_transitions.size(); i++) {
        m_right_transitions.push_back(right_transitions[i]);
    }
}

void regime::add_right_transitions_histogram(const vector< unsigned int > & right_transitions_histogram) {
    if (right_transitions_histogram.size() != m_right_transitions_histogram.size()) {
        cerr << "histogram sizes don't match" << endl;
    }
    for (unsigned int i = 0; i < right_transitions_histogram.size(); i++) {
        m_right_transitions_histogram[i] += right_transitions_histogram[i];
        m_transitions_out += right_transitions_histogram[i];
    }
}

void regime::insert_new_unobserved_regime(const unsigned int & new_unobserved_regime_index) {
    m_right_transitions_histogram.insert(m_right_transitions_histogram.begin() + new_unobserved_regime_index, 0);
    for (unsigned int i = 0; i < m_right_transitions.size(); i++) {
        if (m_right_transitions[i] >= new_unobserved_regime_index && m_right_transitions[i] < numeric_limits< unsigned int >::max()) {
            m_right_transitions[i]++;
        }
    }
}

void regime::remove_unobserved_regime(const unsigned int & regime_to_remove) {
    m_right_transitions_histogram.erase(m_right_transitions_histogram.begin() + regime_to_remove);
    for (unsigned int i = 0; i < m_right_transitions.size(); i++) {
        if (m_right_transitions[i] >= regime_to_remove && m_right_transitions[i] < numeric_limits< unsigned int >::max()) {
            m_right_transitions[i]--;
        }
    }
}

void regime::add_sufficient_statistics(const vector< double > & sufficient_statistics) {
    // check that the sizes of the sufficient statistics match
    if (sufficient_statistics.size() != m_sufficient_statistics.size()) {
        cerr << "sufficient statistics sizes do not match" << endl;
    }
    for (unsigned int i = 0; i < sufficient_statistics.size(); i++) {
        m_sufficient_statistics[i] += sufficient_statistics[i];
    }
}

void regime::remove_sufficient_statistics(const vector< double > & sufficient_statistics) {
    // check that the sizes of the sufficient statistics match
    if (sufficient_statistics.size() != m_sufficient_statistics.size()) {
        cerr << "sufficient statistics sizes do not match" << endl;
    }
    for (unsigned int i = 0; i < sufficient_statistics.size(); i++) {
        m_sufficient_statistics[i] -= sufficient_statistics[i];
    }
}

void regime::add_observations(const unsigned int & trace_index, const double & number_of_observations) {
    m_number_of_observations_in_traces[trace_index] += number_of_observations;
}

void regime::remove_observations(const unsigned int & trace_index, const double & number_of_observations) {
    m_number_of_observations_in_traces[trace_index] -= number_of_observations;
}


void regime::increase_right_indices_greater_than_or_equal_to(const unsigned int & index) {
    // any element of m_right_changepoint_indices that is greater than or equal to index is increased by 1.
    if (m_right_changepoint_indices.size() > 0) {
        unsigned int i = static_cast< unsigned int >(m_right_changepoint_indices.size() - 1);
        bool greater_or_equal = m_right_changepoint_indices[i] >= index;
        bool can_go_lower = true;
        while (greater_or_equal && can_go_lower) {
            m_right_changepoint_indices[i]++;
            i--;
            can_go_lower = i != numeric_limits< unsigned int >::max();
            if (can_go_lower) {
                greater_or_equal = m_right_changepoint_indices[i] >= index;
            }
        }
    }
}

void regime::decrease_right_indices_greater_than(const unsigned int & index) {
    if (m_right_changepoint_indices.size() > 0) {
        unsigned int i = static_cast< unsigned int >(m_right_changepoint_indices.size() - 1);
        bool greater_or_equal = m_right_changepoint_indices[i] >= index;
        bool can_go_lower = true;
        while (greater_or_equal && can_go_lower) {
            m_right_changepoint_indices[i]--;
            i--;
            can_go_lower = i != numeric_limits< unsigned int >::max();
            if (can_go_lower) {
                greater_or_equal = m_right_changepoint_indices[i] >= index;
            }
        }
    }
}

void regime::add_to_similarity_matrix(vector< vector< double > > & similarity_matrix) {
    size_t similarity_matrix_size = similarity_matrix.size();
    for (unsigned int trace_idx_0 = 0; trace_idx_0 < similarity_matrix_size; trace_idx_0++) {
        if (m_number_of_observations_in_traces[trace_idx_0] > 0.5) { // if there are any observations to add for this trace
            for (unsigned int trace_idx_1 = trace_idx_0 + 1; trace_idx_1 < similarity_matrix_size; trace_idx_1++) {
                if (m_number_of_observations_in_traces[trace_idx_1] > 0.5) {
                    similarity_matrix[trace_idx_0][trace_idx_1] += m_number_of_observations_in_traces[trace_idx_0] + m_number_of_observations_in_traces[trace_idx_1];
                    similarity_matrix[trace_idx_1][trace_idx_0] += m_number_of_observations_in_traces[trace_idx_0] + m_number_of_observations_in_traces[trace_idx_1];
                }
            }
        }
    }
}

void regime::add_to_min_proportion_similarity_matrices(vector< vector< double > > & min_proportion_similarity_matrix_0, vector< vector< double > > & min_proportion_similarity_matrix_1, const vector< unsigned long int > & actual_number_of_observations_in_traces) {
    size_t min_proportion_similarity_matrix_size = min_proportion_similarity_matrix_0.size();
    for (unsigned int trace_idx_0 = 0; trace_idx_0 < min_proportion_similarity_matrix_size; trace_idx_0++) {
        if (m_number_of_observations_in_traces[trace_idx_0] > 0.5) { // if there are any observations to add for this trace
            double actual_number_of_observations_0 = static_cast< double >(actual_number_of_observations_in_traces[trace_idx_0]);
            for (unsigned int trace_idx_1 = trace_idx_0 + 1; trace_idx_1 < min_proportion_similarity_matrix_size; trace_idx_1++) {
                if (m_number_of_observations_in_traces[trace_idx_1] > 0.5) { // if trace trace_idx_1 is observed for this regime
                    double actual_number_of_observations_1 = static_cast< double >(actual_number_of_observations_in_traces[trace_idx_1]);
                    min_proportion_similarity_matrix_0[trace_idx_0][trace_idx_1] += m_number_of_observations_in_traces[trace_idx_0] / actual_number_of_observations_0;
                    min_proportion_similarity_matrix_0[trace_idx_1][trace_idx_0] += m_number_of_observations_in_traces[trace_idx_0] / actual_number_of_observations_0;
                    min_proportion_similarity_matrix_1[trace_idx_0][trace_idx_1] += m_number_of_observations_in_traces[trace_idx_1] / actual_number_of_observations_1;
                    min_proportion_similarity_matrix_1[trace_idx_1][trace_idx_0] += m_number_of_observations_in_traces[trace_idx_1] / actual_number_of_observations_1;
                }
            }
        }
    }
}

void regime::add_new_trace() {
    m_number_of_observations_in_traces.push_back(0);
}

void regime::permute_transitions_histogram(const vector< unsigned int > & permutation_vector) {
	vector< unsigned int > histogram_copy = m_right_transitions_histogram;
	for (unsigned int i = 0; i < m_right_transitions_histogram.size(); i++) {
		m_right_transitions_histogram[permutation_vector[i]] = histogram_copy[i];
	}
}

void regime::permute_right_transitions(const vector< unsigned int > & permutation_vector) {
	for (unsigned int i = 0; i < m_right_transitions.size(); i++) {
		if (m_right_transitions[i] < numeric_limits< unsigned int >::max()) {
			m_right_transitions[i] = permutation_vector[m_right_transitions[i]];
		}
	}
}

#endif
