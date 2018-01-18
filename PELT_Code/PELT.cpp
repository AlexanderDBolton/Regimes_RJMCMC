#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <memory>
#include "multiple_processes_regime.h"
using namespace std;

void calculate_F_tau_star_and_tau_1(vector< double > &, vector< unsigned long int > &, double &, unsigned long int &, unsigned long int, mult_process&);
void set_R_tau_star_plus_1(vector< unsigned long int > &, vector< unsigned long int > &, vector< double > &, unsigned long int, const double &, mult_process&);
void write_cp_n_to_file(vector< unsigned long int > &, string, const unsigned long int &);
void print_asterisk_if_needed(unsigned long int, unsigned long int, unsigned long int &, unsigned long int);

int main(int argc, char * argv[]) {
    string data_file = ""; // name of the data file
    string penalty = "BIC"; // can be BIC or AIC
    unsigned long int end_time = 1; // assuming that the process starts at time 0 (with the first transition happening at time 1 for a Markov chain), when does it end?
    unsigned long int diff = 1; // squash multiple observations into a single time
    double beta = 0; // beta will be assigned in mult_process
    double K = 0;
    
    if(argc > 1) {
        ifstream InputParameters(argv[1], ios::in);
        InputParameters >> data_file;
        InputParameters >> penalty;
        InputParameters >> end_time;
        InputParameters >> diff;
    }
    
    std::unique_ptr<mult_process> pm = std::unique_ptr<mult_process>(new mult_process(data_file, penalty, end_time, diff, beta));
    end_time /= diff;
  
    vector< double > F(end_time + 2, 0);
    vector< vector< unsigned long int > > R(end_time + 3);
    vector< vector< unsigned long int > > cp(end_time + 2);
  
    F[0] = -beta;
    R[1].push_back(0);
    unsigned long int number_of_asterisks_already_printed = 0;

    for (unsigned long int tau_star = 1; tau_star <= end_time + 1; tau_star++) {
        unsigned long int tau_1;
        calculate_F_tau_star_and_tau_1(F, R[tau_star], beta, tau_1, tau_star, *pm);
        cp[tau_star] = cp[tau_1];
        cp[tau_star].push_back(tau_1);
        set_R_tau_star_plus_1(R[tau_star + 1], R[tau_star], F, tau_star, K, *pm);
        print_asterisk_if_needed(tau_star, end_time, number_of_asterisks_already_printed, 10);
    }
    
    cout << "\n";
    
    data_file.erase(data_file.end() - 15, data_file.end());
    string cp_n_Filename = data_file + "_cp_n.txt";
    write_cp_n_to_file(cp.back(), cp_n_Filename, diff);
    
    return 0;
}

void calculate_F_tau_star_and_tau_1(vector< double > & F, vector< unsigned long int > & R_tau_star, double & beta, unsigned long int & tau_1, unsigned long int tau_star, mult_process & pm) {
    tau_1 = R_tau_star[0];
    double min = F[tau_1] + pm.cost(tau_1, tau_star) + beta;
    for (unsigned int index = 1; index < R_tau_star.size(); index++) {
        unsigned long int tau = R_tau_star[index];
        double cost = F[tau] + pm.cost(tau, tau_star) + beta;
        if (cost < min) {
            min = cost;
            tau_1 = tau;
        }
    }
    F[tau_star] = min;
}

void set_R_tau_star_plus_1(vector< unsigned long int > & R_tau_star_plus_1, vector< unsigned long int > & R_tau_star, vector< double > & F, unsigned long int tau_star, const double & K, mult_process & pm) {
    for (unsigned int index = 0; index < R_tau_star.size(); index++) {
        if (F[R_tau_star[index]] + pm.cost(R_tau_star[index], tau_star) + K <= F[tau_star]) {
            R_tau_star_plus_1.push_back(R_tau_star[index]);
        }
    }
    R_tau_star_plus_1.push_back(tau_star);
}

void write_cp_n_to_file(vector< unsigned long int > & cp_n, string output_filename, const unsigned long int & diff) {
    ofstream OutputStream(output_filename, ios::out);
    OutputStream << setiosflags(ios::fixed);
    for (unsigned long int i = 1; i < cp_n.size(); i++) {
        OutputStream << cp_n[i] * diff << endl;
    }
    OutputStream.close();
}

void print_asterisk_if_needed(unsigned long int x, unsigned long int n, unsigned long int & number_already_printed, unsigned long int num_asterisks) {
    if (x > n * number_already_printed / num_asterisks) {
        cout << "*";
        cout.flush();
        number_already_printed++;
    }
}
